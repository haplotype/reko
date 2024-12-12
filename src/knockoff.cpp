#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/Jacobi>
#include <vector>
#include <algorithm>
#include <random>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

#include <unistd.h>
#include "thread_pool.h" // Include the thread pool header

// Function to center and scale each column of a matrix (similar to R's scale function)
Eigen::MatrixXd scale(const Eigen::MatrixXd& X) {
    Eigen::MatrixXd X_scaled = X;
    for (int j = 0; j < X.cols(); ++j) {
        double mean = X.col(j).mean();
        double stddev = std::sqrt((X.col(j).array() - mean).square().sum() / (X.rows() - 1));
	if(stddev > 1e-10) 
           X_scaled.col(j) = (X.col(j).array() - mean) / stddev;
	else
           X_scaled.col(j).setZero(); 
    }
    return X_scaled;
}

// Function to compute the covariance matrix between two matrices
Eigen::MatrixXd covariance(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y) {
    Eigen::MatrixXd X_centered = X.rowwise() - X.colwise().mean();
    Eigen::MatrixXd Y_centered = Y.rowwise() - Y.colwise().mean();
    return (X_centered.transpose() * Y_centered) / static_cast<double>(X.rows() - 1);
}

// The sphere_vec function with Step 4 parallelized using a thread pool
std::vector<Eigen::MatrixXd> sphere_vec(const Eigen::MatrixXd& xin, int num = 1, int nth = 8) {
    int ni = xin.rows();
    int np = xin.cols();

    // Step 1: Append a diagonal matrix to xin
    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Identity(np, np) * std::sqrt(static_cast<double>(np));
    Eigen::MatrixXd new_xin(ni + np, np);
    new_xin << scale(xin), diagMatrix;

    // Step 2: Scale the combined matrix
    Eigen::MatrixXd mX = scale(new_xin);

    // Step 3: Perform QR decomposition
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(mX);
    Eigen::MatrixXd xQ = qr.householderQ();
    Eigen::MatrixXd xR = qr.matrixQR().triangularView<Eigen::Upper>();

    Eigen::MatrixXd tR = xR;
    //note xR and tR is (ni+np) by np, not np by np;
//    std::cout << xR.rows() << " xR " << xR.cols() << std::endl; 
//    std::cout << xR.topRows(20) << std::endl; 

    // Step 4: Modify tR using projections (Parallelized using thread pool)
    int num_cols = tR.cols();

    // Create a thread pool with a specified number of threads
//    size_t numThreads = std::thread::hardware_concurrency(); // Use the number of hardware threads
    long numThreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (numThreads < 1) 
    {
        numThreads = nth; // Fallback to nth (8) if detection fails
    }
    ThreadPool pool(numThreads);

    // Since Eigen matrices are not thread-safe for writes, we ensure no overlapping writes
    for (int j = 0; j < num_cols; ++j) {
        pool.enqueue([&, j]() {
            // Exclude column j from tR
	    Eigen::VectorXd x = xR.col(j).topRows(np);
            Eigen::MatrixXd temp(np, np-1); 
	    int c = 0; 
	    for(int k = 0; k < np; k++)
	    {
	        if(k == j) continue; 
                temp.col(c) = xR.col(k).topRows(np); 
		c++; 
//		std::cout << c << std::endl;
	    }
//	    if (j > 0) {
//               temp.leftCols(j) = xR.leftCols(j);
//	    }
//	    if(j < np-1) {
//	       temp.rightCols(np-1-j) = xR.rightCols(np-1-j);
//	    }
//	    temp.col(np-1).setZero(); 
//	    temp(np-1,np-1)=1;
//	    std::cout << temp << std::endl; 
            Eigen::MatrixXd temp_copy = temp; 

	    for(int k = 0; k < np-1; k++) {
                if(fabs(temp(k+1,k)) > 1e-10) {
		    double a = temp(k,k); 
		    double b = temp(k+1,k); 
		    double r = sqrt(a*a+b*b); 
		    double cos = a / r; 
		    double sin = b / r; 

//		    // Compute Givens rotation coefficients
//		    Eigen::JacobiRotation<double> G;
//		    G.makeGivens(a,b);
		    for(int c = 0; c < np-1; c++) {
			double temp_a = cos * temp(k,c) + sin * temp(k+1,c);
			double temp_b = -sin * temp(k,c) + cos * temp(k+1,c);
			temp(k,c) = temp_a; 
			temp(k+1,c) = temp_b; 
		    }

		    // Apply rotation to elements i-1 and i of x
//                    G.apply(x(k), x(k+1));
		    double temp_x = cos * x(k) + sin * x(k+1);
                    x(k+1) = -sin * x(k) + cos * x(k+1);
                    x(k) = temp_x;
		}
//	        std::cout << temp << std::endl; 
	    }

//	    std::cout << temp << std::endl; 
//	    std::cout << x << std::endl; 

	    // Solve the system temp * beta = x
            Eigen::MatrixXd temp2=temp.topRows(np-1); 
            Eigen::VectorXd x2 = x.topRows(np-1); 
            Eigen::VectorXd beta = temp2.triangularView<Eigen::Upper>().solve(x2);
//	    beta(np-1) = 0; 
//	    std::cout << beta << std::endl; 
            // Compute fj = temp * beta
            Eigen::VectorXd fj = temp_copy * beta;
//	    std:: cout << fj << std::endl; 

//	    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_temp(temp);
//            Eigen::VectorXd beta = qr_temp.solve(xR.col(j));
//	    std::cout << "beta.size " << beta.size() << std::endl; 
//            Eigen::VectorXd fj = temp * beta;
//	    std:: cout << "fj.size " << fj.size() << std::endl; 

            // Update tR
            tR.col(j).topRows(np) = fj + fj - xR.col(j).topRows(np);
//	    std:: cout << tR.col(j).topRows(np) << std::endl; 
//            tR.col(j) = fj + fj - xR.col(j);
        });
    }

    // Wait for all tasks to complete
    pool.wait();

    // Step 5: Compute mY
    Eigen::MatrixXd mY = xQ * tR;
//    std:: cout << mY.topRows(np*2) << std::endl; 
   

    // Step 6: Compute mW
    Eigen::MatrixXd mW = (mX - mY) / 2.0;

    // Step 7: Compute covariance matrices
    Eigen::MatrixXd mT = covariance(mW, mX);
    Eigen::MatrixXd mW_scaled = scale(mW);
    Eigen::MatrixXd tdt = covariance(mW_scaled, mW_scaled);

    // Step 8: Eigen decomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(tdt);
    Eigen::VectorXd lam = eigensolver.eigenvalues();
    Eigen::MatrixXd Q = eigensolver.eigenvectors();

    // Step 9: Compute alpha
    double max_lam = lam.maxCoeff();
    double alpha = 1.0 / max_lam;
    alpha = std::min(1.0, alpha);

    // Step 10: Compute B
    Eigen::VectorXd B = 4 * alpha * (1 - alpha * lam.array());
    B = B.cwiseMax(0.0).array().sqrt();

    // Step 11: Prepare the output container
    std::vector<Eigen::MatrixXd> ko(num, Eigen::MatrixXd(ni, np));

    // Step 12: Compute the mean matrix
    Eigen::MatrixXd mean = (1.0 - alpha) * mX.topRows(ni) + alpha * mY.topRows(ni);

    // Step 13: Random number generation setup
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> normal(0.0, 1.0);

    // Compute the sqrt of the diagonal of mT
    Eigen::VectorXd sqrt_mT_diag = mT.diagonal().array().sqrt();
    // Step 14: Generate random matrices and compute ko
    for (int j = 0; j < num; ++j) {
        // Generate a random matrix with standard normal entries
        Eigen::MatrixXd uu = Eigen::MatrixXd::NullaryExpr(ni, np, [&]() { return normal(generator); });
        uu = scale(uu);

        // Compute ko[j]
        ko[j] = mean + uu * B.asDiagonal() * Q.transpose() * sqrt_mT_diag.asDiagonal();
        ko[j] = scale(ko[j]);
    }

    return ko;
}

// Function to read a numerical matrix from a file
Eigen::MatrixXd readMatrixFromFile(const std::string& filename, int& n, int& p) {
    std::ifstream infile(filename);
    std::string line;
    std::vector<double> values;
    int cols = -1;
    n = 0;

    if (!infile) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    while (std::getline(infile, line)) {
        // Ignore empty lines
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::vector<double> rowValues;
        double val;
        while (iss >> val) {
            rowValues.push_back(val);
        }

        // Skip lines with no data
        if (rowValues.empty()) continue;

        // Determine the number of columns
        if (cols == -1) {
            cols = static_cast<int>(rowValues.size());
        } else if (static_cast<int>(rowValues.size()) != cols) {
            throw std::runtime_error("Inconsistent number of columns in row " + std::to_string(n + 1));
        }

        values.insert(values.end(), rowValues.begin(), rowValues.end());
        ++n;
    }

    p = cols;

    // Check if the file was empty or had inconsistent data
    if (n == 0 || p == 0) {
        throw std::runtime_error("No data found in file or file is incorrectly formatted.");
    }

    // Create an Eigen matrix and fill it with the values
    Eigen::MatrixXd mat(n, p);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            mat(i, j) = values[i * p + j];
        }
    }

    return mat;
}

// Main function to handle command-line arguments and perform the computations
int main(int argc, char* argv[]) {
    std::string filename;
    int num = 1;
    int nth = 8; 
    std::string output_prefix("out");

    // Simple command-line argument parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-f" && i + 1 < argc) {
            filename = argv[++i];
        } else if (arg == "-t" && i + 1 < argc) {
            nth = std::stoi(argv[++i]);
        } else if (arg == "-c" && i + 1 < argc) {
            num = std::stoi(argv[++i]);
        } else if (arg == "-o" && i + 1 < argc) {
            output_prefix = argv[++i];
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            std::cerr << "Usage: " << argv[0] << " -f filename -c number -t thread -o output_prefix\n";
            return 1;
        }
    }

    if (filename.empty() || output_prefix.empty()) {
        std::cerr << "Filename (-f) and output prefix (-o) must be specified.\n";
        std::cerr << "Usage: " << argv[0] << " -f filename -c number -t thread -o output_prefix\n";
        return 1;
    }

    // Read the matrix from the file
    int n = 0, p = 0;
    Eigen::MatrixXd xin;
    try {
        xin = readMatrixFromFile(filename, n, p);
    } catch (const std::exception& e) {
        std::cerr << "Error reading matrix from file: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "number of samples (rows) = " << n << std::endl; 
    std::cout << "number of features (cols) = " << p << std::endl; 

    // Generate the knockoff matrices
    std::vector<Eigen::MatrixXd> ko_matrices;
    try {
        ko_matrices = sphere_vec(xin, num, nth);
    } catch (const std::exception& e) {
        std::cerr << "Error generating knockoff matrices: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "copies of knockoffs = " << num << std::endl; 
    // Write each matrix to a separate file
    for (int i = 0; i < num; ++i) {
        std::string output_filename = output_prefix + "." + std::to_string(i + 1) + ".ko";
        std::ofstream outfile(output_filename);
        if (!outfile) {
            std::cerr << "Error opening file for writing: " << output_filename << std::endl;
            return 1;
        }
        outfile << ko_matrices[i] << "\n";
        outfile.close();
    }

    return 0;
}


#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/Jacobi>
#include <vector>
#include <algorithm>
#include <random>
#include <omp.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <functional>

#include <cstdlib>
#include <unistd.h>
#include "thread_pool.h" // Include the thread pool header

#define VERSION "0.314" //17 MARCH 2025

void rSVD2G(const Eigen::MatrixXd &G, int k, int p, int q, 
                  Eigen::MatrixXd &U_A, Eigen::VectorXd &evals_A, int randseed)
{
    int n = (int)G.rows();
    int m = (int)G.cols();
    int l = k + p;
    assert(m>l); 

//    srand(time(NULL)); 
//    int randseed = rand() % 1000001; 
    std::mt19937 gen(randseed);
    std::normal_distribution<> dist(0.0,1.0);

    Eigen::MatrixXd Omega(n, l);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            Omega(i,j) = dist(gen);
        }
    }

    // (1/m GG^t - I);   
    Eigen::MatrixXd W = G.transpose() * Omega; 
    Eigen::MatrixXd Y_A = G * W;               

    for (int i = 0; i < q; i++) {
        Eigen::MatrixXd W_A = G.transpose() * Y_A;
        Y_A = G * W_A;
    }

    Eigen::HouseholderQR<Eigen::MatrixXd> qrA(Y_A);
    Eigen::MatrixXd Q_A = qrA.householderQ() * Eigen::MatrixXd::Identity(n, l);

    Eigen::MatrixXd B_A = Q_A.transpose() * G;
    Eigen::JacobiSVD<Eigen::MatrixXd> svdA(B_A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd U_B_A = svdA.matrixU();
    Eigen::VectorXd S = svdA.singularValues();

    U_A = Q_A * U_B_A; 
    U_A.conservativeResize(n, k);

    evals_A = S.head(k).array().square();
}
//
// Function to solve the secular equation for updated eigenvalues using bisection method
double solveSecularEquation(const Eigen::VectorXd& old_eigenvalues, const Eigen::VectorXd& v, int j) {
    double lower = (j == old_eigenvalues.size()-1) ? 0 : old_eigenvalues(j + 1);
    double upper = old_eigenvalues(j);
    double tolerance = 1e-8;
    double lambda = 0.0;

    while (upper - lower > tolerance) {
        lambda = (upper + lower) / 2.0;
        double sum = 0.0;
        for (int k = 0; k < old_eigenvalues.size(); ++k) {
            sum += v(k) * v(k) / (old_eigenvalues(k) - lambda);
        }
        if (1 - sum > 0) {
            lower = lambda;
        } else {
            upper = lambda;
        }
    }
    return lambda;
}

// function to update the eigenvectors
void updateeigenvectors(const Eigen::MatrixXd& uxx, const Eigen::VectorXd& dxx, const Eigen::VectorXd& xj, Eigen::MatrixXd& uj) {
    Eigen::VectorXd v = uxx.transpose() * xj; 
    for (int i = 0; i < uj.cols(); ++i) {
//	std::cout << "update eigenvector " << i << std::endl; 
//        uj.col(i) = uxx.col(i) + v(i) * (uxx * (eigen::matrixxd::identity(n, n).col(i) / (dxx(i) - dj(i))));
	double dji = solveSecularEquation(dxx, v, i);
	Eigen::VectorXd diff = dxx.array() - dji;
	Eigen::VectorXd dinv = diff.array().inverse(); 
	Eigen::VectorXd w = dinv.array() * v.array(); 
	Eigen::VectorXd u = uxx * w; 
	uj.col(i) = u / u.norm(); 
    }

}

void updateeigenvectors2(const Eigen::MatrixXd& uxx, const Eigen::VectorXd& dxx, const Eigen::VectorXd& dj,const Eigen::VectorXd& xj, Eigen::MatrixXd& uj) {
    Eigen::VectorXd v = uxx.transpose() * xj; 
    for (int i = 0; i < uj.cols(); ++i) {
//	std::cout << "update eigenvector " << i << std::endl; 
//        uj.col(i) = uxx.col(i) + v(i) * (uxx * (eigen::matrixxd::identity(n, n).col(i) / (dxx(i) - dj(i))));
	double dji = dj(i);
	Eigen::VectorXd diff = dxx.array() - dji;
	Eigen::VectorXd dinv = diff.array().inverse(); 
	Eigen::VectorXd w = dinv.array() * v.array(); 
	Eigen::VectorXd u = uxx * w; 
	uj.col(i) = u / u.norm(); 
    }

}


//int update_lead_eigenvector(const Eigen::MatrixXd & X, const Eigen::VectorXd & x, Eigen::VectorXd & v1) {
//    int steps = 0; 
//    while(steps < 100) {
//	Eigen::VectorXd u = X * (X.transpose() * v1); 
//	Eigen::VectorXd ux = x * x.dot(v1); 
//	u = u - ux;
//	double unorm = u.norm(); 
//	u = u / unorm; 
//	double dnorm = (u-v1).norm(); 
////	    std::cout << v1.norm() << " " << unorm << " " << dnorm <<std::endl; 
//	if(dnorm < 1e-5) break; 
//	v1 = u; 
//	steps++; 
//    }
////    std::cout << std::endl;
//    return steps; 
//}

double df(double x, const Eigen::VectorXd & yQ, const Eigen::VectorXd & D) {

//    w=which(D>1e-6);
//    y2=yQ[w];
//    D2=D[w];
//    H=D2*eta+1;
//    H2=H^2;
//    return(sum((H-y2^2)*D2/H2));  
    double sum = 0; 
    for(unsigned int i = 0; i < D.size(); i++)
    {
	if(D(i) < 1e-6) continue; 
	double y2=yQ(i); 
	double D2=D(i); 
        double H = x * D2 + 1; 
        double H2 = H * H; 
        sum += (H - y2 * y2) * D2 / H2; 
    }

    return sum;                                                                 
}

double bisect(std::function<double(double)> func, double a, double b, double tol)
{
    double res=1;  
    if(func(a) * func(b) < 0) {
	double c = a; 
	while ((b-a) > tol) {
	    c = (a+b) / 2.0; 
	  if(std::abs(func(c)) < tol) {
	      break; 
	  }             

	  if(func(c) * func(a) < 0) {
	      b=c; 
	  }
	  else {
	      a=c; 
	  }
	}
	res = c; 
    }
    return std::sqrt(res); 
}

// Function to perform SVD and return U and D^2
void computeSVD(const Eigen::MatrixXd& X, Eigen::MatrixXd& U, Eigen::VectorXd& D2) {
    // Perform SVD: X = U * S * V^T
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    U = svd.matrixU();
//    V = svd.matrixV();
    Eigen::VectorXd S = svd.singularValues();
    
    // Compute D^2 (eigenvalues of XX^T)
    D2 = S.array().square();
}


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
int sigma_eb(const Eigen::MatrixXd& xin, Eigen::VectorXd &vsig, int nth, int randseed, int ebpc) {
//    int ni = xin.rows();
    int np = xin.cols();

//    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Identity(np, np) * std::sqrt(1);
//    Eigen::MatrixXd new_xin(ni + np, np);
//    new_xin << scale(xin), diagMatrix;

    Eigen::MatrixXd mX = scale(xin);
    
    ///////////////////////////////////////////////////////
    //Empirical Bayes to estimate optimal sigma for each Householder reflection. 
    //EB method is base on eigen decomposition of X[,-j]X[,-j]^t. 
    // Perform SVD to get U and D^2 for XX^T
    Eigen::MatrixXd Uxx;
    Eigen::VectorXd Dxx;
    if(ebpc == 0) 
       computeSVD(mX, Uxx, Dxx);
    else
//    if(np > ebpc + 10) 
       rSVD2G(mX, ebpc, 10, 2, Uxx, Dxx, randseed); 

    long numThreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (numThreads < 1) 
    {
        numThreads = nth; // Fallback to nth (8) if detection fails
    }
    ThreadPool pool(numThreads);

    for (int j = 0; j < mX.cols(); j++) {
        pool.enqueue([&, j]() {
	// reverse rank-one update to obtain eigenvalues for X_[,-j]X[,-j]^t; 
	    Eigen::VectorXd Xj = mX.col(j);
	    Eigen::VectorXd Dj(Uxx.cols());
	    for (int k = 0; k < Uxx.cols(); ++k) {
		Eigen::VectorXd v = Uxx.transpose() * Xj; 
		Dj(k) = solveSecularEquation(Dxx, v, k);
    //	    std::cout << k << "\t " << Dj(k) << std::endl; 
	    }

	    // Update the eigenvectors
	    Eigen::MatrixXd Uj = Uxx;
//   	    updateeigenvectors2(Uxx, Dxx, Dj, Xj, Uj);

	    Eigen::VectorXd yQ = Uj.transpose() * Xj; 
	    auto func = [&yQ, &Dj](double x) {
		return df(x,yQ,Dj);
	    }; 
	    double inv_sig = 1.0 / bisect(func, 0, 1, 1e-6); 
            vsig(j) = inv_sig; 

        });
    }
    // Wait for all tasks to complete
    pool.wait();

    return 1;
}


std::vector<Eigen::MatrixXd> sphere_inv(const Eigen::MatrixXd& xin, const Eigen::VectorXd & vsig, double lambda, int ebpc, int npc, int randseed, int ncopy, int nth) {
    int ni = xin.rows();
    int np = xin.cols();
//    if(npc > np) npc = np; 
//    if(npc == 0) 
//	npc = (int) (np * 0.2); 
                      
//    double inv_sig = 1; 
//    if(rint(inv_sigma) == 0) inv_sig = std::sqrt(np); 
//    else if(rint(inv_sigma) == 1) inv_sig = std::sqrt(ni); 
//    else inv_sig = inv_sigma; 

//    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Identity(np, np) * inv_sig;
    Eigen::MatrixXd diagMatrix = vsig.asDiagonal(); 
    Eigen::MatrixXd new_xin(ni + np, np);
    new_xin << scale(xin), diagMatrix;

    Eigen::MatrixXd mX = scale(new_xin);
    Eigen::MatrixXd mY = mX;

    Eigen::MatrixXd U;
    Eigen::VectorXd D2;
    if(ebpc == 0) 
       computeSVD(mX, U, D2);
    else //    if(np > npc + 10) 
       rSVD2G(mX, ebpc, 10, 2, U, D2, randseed); 

//    for(int j = 0; j < D2.size(); j++)
//	std::cout << D2(j) << std::endl; 

//    int wh = 0; 
//    Eigen::VectorXd cumsum = D2;
//    for(int j = 1; j < D2.size(); j++)
//    {
//	cumsum(j) = cumsum(j-1) + D2(j);  
//    }
//    for(int j = 0; j < D2.size(); j++)
//    {
//	cumsum(j) /= cumsum(D2.size()-1);  
//	if(cumsum(j) < lambda) 
//	    wh = j;
//    }
//    std::cout << wh << std::endl; 
//    for(int j = 0; j < cumsum.size(); j++)
//       std::cout << cumsum(j) << " "; 
//    std::cout << std::endl; 
//
//    npc = wh; 

    Eigen::MatrixXd X = U.leftCols(npc);
//    std::cout << X.cols() << " " << X.rows() << std::endl;

    long numThreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (numThreads < 1) 
    {
        numThreads = nth; // Fallback to nth (8) if detection fails
    }

    ThreadPool pool(numThreads);
    for(int j = 0; j < np; j++)
    {
        pool.enqueue([&, j]() {
//	    Eigen::MatrixXd X1(ni+np,np-1); 
//	    X1 << mX.leftCols(j), mX.rightCols(np-j-1); 
//	    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mX, Eigen::ComputeThinU | Eigen::ComputeThinV);
////	    Eigen::VectorXd S = svd.singularValues();
//	    Eigen::MatrixXd U = svd.matrixU();
//            Eigen::MatrixXd X = U.leftCols(npc);

//	    Eigen::VectorXd bhat = X.transpose() * mX.col(j); 
////	    std::cout << X.norm() << " " << mX.col(j).norm() << " " << bhat << std::endl; 
//	    Eigen::VectorXd fj = X * bhat; 
//	    mY.col(j) = fj + fj - mX.col(j); 

	    Eigen::VectorXd y = mX.col(j); 
            Eigen::MatrixXd Uj = X;
//            updateeigenvectors(U, D2, y, Uj);  
	    Eigen::VectorXd bhat = Uj.transpose() * y; 
	    Eigen::VectorXd fj = Uj * bhat;   
//	    std::cout << v1.norm() << " " << y.norm() << " " << bhat << " " <<  fj.norm() << std::endl; 
	    mY.col(j) = fj + fj - y; 

        });
    }
    // Wait for all tasks to complete
    pool.wait();

   
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
    std::vector<Eigen::MatrixXd> ko(ncopy, Eigen::MatrixXd(ni, np));

    // Step 12: Compute the mean matrix
    Eigen::MatrixXd mean = (1.0 - alpha) * mX.topRows(ni) + alpha * mY.topRows(ni);

    // Step 13: Random number generation setup
//    std::random_device rd;
    std::mt19937 generator(randseed+1);
    std::normal_distribution<double> normal(0.0, 1.0);

    // Compute the sqrt of the diagonal of mT
    Eigen::VectorXd sqrt_mT_diag = mT.diagonal().array().sqrt();
    // Step 14: Generate random matrices and compute ko
    for (int j = 0; j < ncopy; ++j) {
        // Generate a random matrix with standard normal entries
        Eigen::MatrixXd uu = Eigen::MatrixXd::NullaryExpr(ni, np, [&]() { return normal(generator); });
        uu = scale(uu);

        // Compute ko[j]
        ko[j] = mean + uu * B.asDiagonal() * Q.transpose() * sqrt_mT_diag.asDiagonal();
        ko[j] = scale(ko[j]);
    }

    return ko;
}


// The sphere_vec function with Step 4 parallelized using a thread pool
std::vector<Eigen::MatrixXd> sphere_vec(const Eigen::MatrixXd& xin, const Eigen::VectorXd & vsig, int randseed, int ncopy, int nth) {
    int ni = xin.rows();
    int np = xin.cols();

    Eigen::MatrixXd diagMatrix = vsig.asDiagonal(); 
    Eigen::MatrixXd new_xin(ni + np, np);
    new_xin << scale(xin), diagMatrix;

//    Eigen::MatrixXd diagMatrix = Eigen::MatrixXd::Identity(np, np) * inv_sigma;
//    Eigen::MatrixXd new_xin(ni + np, np);
//    new_xin << scale(xin), diagMatrix;

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
    std::vector<Eigen::MatrixXd> ko(ncopy, Eigen::MatrixXd(ni, np));

    // Step 12: Compute the mean matrix
    Eigen::MatrixXd mean = (1.0 - alpha) * mX.topRows(ni) + alpha * mY.topRows(ni);

    // Step 13: Random number generation setup
//    std::random_device rd;
//    std::mt19937 generator(rd());
    std::mt19937 generator(randseed+1);
    std::normal_distribution<double> normal(0.0, 1.0);

    // Compute the sqrt of the diagonal of mT
    Eigen::VectorXd sqrt_mT_diag = mT.diagonal().array().sqrt();
    // Step 14: Generate random matrices and compute ko
    for (int j = 0; j < ncopy; ++j) {
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

int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Version: %s\n", VERSION);
	fprintf(stderr, "Usage:   reko -f filename -c number -t thread -o output_prefix [-b:c:f:o:r:s:t:]\n");
        fprintf(stderr, "Options: \n");
        fprintf(stderr, "         -b int        number of eigenpairs for rSVD in approximate EB estimates [0]\n");
        fprintf(stderr, "         -c int        copies of knockoffs to be generated [1]\n");
        fprintf(stderr, "         -f str        input file contain features nxp\n");
//	fprintf(stderr, "         -l flt        used to select leading PCs for reflection [0.7]\n");
//	fprintf(stderr, "         -n int        number of top PCs to approximate reflection [0]\n");
	fprintf(stderr, "         -o str        output prefix [out]\n");
	fprintf(stderr, "         -r int        random seedL\n");                            
//	fprintf(stderr, "         -p int        number of top PCs to approximate reflection [0]\n");
	fprintf(stderr, "         -s char       e: empirical Baeys estimates; p: inv_sigma = sqrt(p) [e]\n");
	fprintf(stderr, "            flt        inv_sigma where b~MNV(0,V) and V = inv_sigma I_p \n");
//	fprintf(stderr, "         -u flt        threshold to pick columns of U in SVD \n");
	fprintf(stderr, "         -t int        number of threads [8]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Bug report: Yongtao Guan <ytguan@gmail.com>\n\n");
	return 1;
}

// Main function to handle command-line arguments and perform the computations
int main(int argc, char* argv[]) {
    std::string filename;
    int ncopy = 10;
    int nth = 8; 
    char s_char = 'e'; 
    //s_char can be 'e', 'p','n'; 
    double s_val = 0;    //inv_sigma; 
    int npc = 0; 
    int ebpc = 0; 
    std::string output_prefix("out");

    double lambda = 0.5; 
    int randseed = -1; 

    // Simple command-line argument parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-f" && i + 1 < argc) {
            filename = argv[++i];
        } else if (arg == "-t" && i + 1 < argc) {
            nth = std::stoi(argv[++i]);
        } else if (arg == "-b" && i + 1 < argc) {
            ebpc = std::stoi(argv[++i]);
        } else if (arg == "-c" && i + 1 < argc) {
            ncopy = std::stoi(argv[++i]);
        } else if (arg == "-s" && i + 1 < argc) {
	    std::string next_arg = argv[i+1]; 
	    char* end; 
	    double val = strtod(next_arg.c_str(), &end); 
	    if(end != next_arg.c_str() && *end == 'e') {
                s_val = val;
		i++; 
	    } else if (next_arg.length() == 1) {
		s_char = next_arg[0]; 
		i++;
	    }
        } else if (arg == "-r" && i + 1 < argc) {
            randseed = std::stoi(argv[++i]);
//        } else if (arg == "-n" && i + 1 < argc) {
//            npc = std::stoi(argv[++i]);
        } else if (arg == "-l" && i + 1 < argc) {
            lambda = std::stod(argv[++i]);
        } else if (arg == "-o" && i + 1 < argc) {
            output_prefix = argv[++i];
        } else if (arg == "-p" && i + 1 < argc) {
            npc = std::stoi(argv[++i]);
        } else if (arg == "-h" && i + 1 == argc) {
	    usage(); 
	    return 1; 
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
	    usage(); 
            return 1;
        }
    }

    if (filename.empty() || output_prefix.empty()) {
        std::cerr << "Filename (-f) and output prefix (-o) must be specified.\n";
	usage(); 
        return 1;
    }

    srand(time(NULL)); 
    if(randseed < 0) randseed = rand() % 1000001; 

    std::cout << "random seed: " << randseed << std::endl; 
    if(npc == 0) 
	std::cout << "reflection over all other features." << std::endl; 
    else 
	std::cout << "reflection over all top PCs of all features" << std::endl; 
    std::cout << "number of top PCs for EB estimates: " << ebpc << std::endl; 
    std::cout << "number of top PCs for reflection: " << npc << std::endl; 
    std::cout << "number of knockoff copies: " << ncopy << std::endl; 
    std::cout << s_char << " : " << s_val << std::endl; 
    if(s_char == 'e') 
	std::cout << "data augmentation by priors estimated with empirical Bayes method." << std::endl; 
    else if(s_char == 'p') 
	std::cout << "data augmentation by diag matrix of sqrt(p)." << std::endl; 
    else if(s_char == 'n') 
	std::cout << "data augmentation by diag matrix of sqrt(n)." << std::endl; 
    else if(s_char == 'o') 
	std::cout << "data augmentation by diag matrix of sqrt(sqrt(np))." << std::endl; 
    else 
	std::cout << "data augmentation by diag matrix of " << s_val << "." << std::endl; 


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


    ///////////////////////////////////////////////////////
    // Generate the knockoff matrices
    std::vector<Eigen::MatrixXd> ko_matrices;
    try {
	  Eigen::VectorXd vsig(xin.cols()); 
	  if(s_char == 'e')
	  {
	      sigma_eb(xin, vsig, nth, randseed, ebpc);
	      for(int j = 0; j < vsig.size(); j++)
		  std::cout << vsig(j) << std::endl; 
//	      double sum = 0; 
//	      for(int j = 0; j < vsig.size(); j++)
//		  sum += 1.0/vsig(j); 
//	      sum = vsig.size() / sum; 
//	      for(int j = 0; j < vsig.size(); j++)
//		  vsig(j) = sum; 
//	      std::cout << sum << std::endl; 
	  }
	  else {
	      if(s_char == 'p') 
		  s_val = std::sqrt(xin.cols()); 
	      else if(s_char == 'n') 
		  s_val = std::sqrt(xin.rows()); 
	      else if(s_char == 'o') 
		  s_val = std::sqrt(std::sqrt(xin.rows() * xin.cols())); 

              for(int j = 0; j < vsig.size(); j++)
		  vsig(j) = s_val;  
	  }
	  if(npc > 0) 
              ko_matrices = sphere_inv(xin, vsig, lambda, ebpc, npc, randseed, ncopy, nth);
	  else 
              ko_matrices = sphere_vec(xin, vsig, randseed, ncopy, nth);
    } catch (const std::exception& e) {
        std::cerr << "Error generating knockoff matrices: " << e.what() << std::endl;
        return 1;
    }

    // Write each matrix to a separate file
    for (int i = 0; i < ncopy; ++i) {
        std::string output_filename = output_prefix + "." + std::to_string(i + 1) + ".ko";
//        std::ofstream outfile(output_filename);
//        if (!outfile) {
//            std::cerr << "Error opening file for writing: " << output_filename << std::endl;
//            return 1;
//        }
//        outfile << ko_matrices[i] << std::endl;
//        outfile.close();

	FILE * outfile = fopen(output_filename.c_str(),"w"); 
	if(!outfile) {
	    std::cerr << "Error opening file" << std::endl; 
	    return 1; 
	}
	for (int r=0; r<ko_matrices[i].rows(); ++r) {
	    for (int c=0; c<ko_matrices[i].cols(); ++c) {
		fprintf(outfile, "%8.7f" , ko_matrices[i](r,c)); 
		if(c < ko_matrices[i].cols()-1) {
		    fprintf(outfile, " "); 
		}
	    }
	    fprintf(outfile, "\n"); 
	}
	fclose(outfile); 
    }

    return 0;
}


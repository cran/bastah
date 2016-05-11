bastah.scalreg <- function(X, y){
   nX = dim(X)[1];
   pX = dim(X)[2];
   lam0 = sqrt(2 * log(pX) / nX);

   if (requireNamespace("rPython", quietly = TRUE)) {
      numpy = tryCatch(rPython::python.exec("import numpy"), error = function(e)("Error"));
      sklearn = tryCatch(rPython::python.exec("from sklearn import linear_model"), error = function(e)("Error"));
   } else {
      numpy = "Error";
      sklearn = "Error";
   }

   if(is.null(numpy) & is.null(sklearn)){
      rPython::python.assign("x", X);
      rPython::python.assign("y", y);
      rPython::python.exec("x = numpy.asarray(x)");
      rPython::python.exec("y = numpy.asarray(y)");
      rPython::python.exec("itr = 8 * x.shape[1]");
      rPython::python.exec("alphas, active, coef = linear_model.lars_path(x, y, method='lasso', max_iter=itr)");

      lambda = rPython::python.get("alphas.tolist()");
      lambda = lambda * nX;

      beta = rPython::python.get("coef.tolist()");
      beta <- matrix(unlist(beta), nrow = length(lambda), byrow = FALSE);
      lambda = lambda[1:(length(lambda) - 1)];

      objlasso <- list(call = match.call(), type = "LASSO", df=NULL, lambda=lambda,R2 = 1, RSS = NULL, Cp = NULL,
                       actions = NULL, entry = NULL, Gamrat = NULL, arc.length = NULL, Gram = NULL,
                       beta = beta, mu = 0, normx = rep(1, pX), meanx = rep(0, pX));
      rm(beta, lambda);
      class(objlasso) <- "lars"
   }else{
      if(requireNamespace("rPython", quietly = TRUE)){
         warning("installing numpy and sklearn can speed up the computation");
      }else{
         warning("package rPython can speed up the computation using numpy and sklearn");
      }
      objlasso = lars(X, y, type = "lasso", intercept = FALSE, normalize = FALSE, use.Gram = FALSE);
   }
   rm(numpy, sklearn, pX);
   #=================================================================
   #Code taken from scalereg package
   sigmaint = 0.1;
   sigmanew = 5;
   flag = 0;
   while(abs(sigmaint-sigmanew) > 0.0001 & flag <= 100){
      flag = flag + 1;
      sigmaint = sigmanew;
      lam = lam0 * sigmaint;
      hy = predict.lars(objlasso, X, s = lam * nX, type = "fit", mode = "lambda")$fit;
      sigmanew = sqrt(mean((y - hy) ^ 2));
   }
   rm(flag, lam0, sigmaint);

   hbeta = predict.lars(objlasso, X, s = lam * nX, type = "coefficients", mode = "lambda")$coef;
   hy = predict.lars(objlasso, X, s = lam * nX, type = "fit", mode = "lambda")$fit;
   est = list(hsigma = sigmanew, coefficients = hbeta, residuals = y - hy, fitted.values = hy);
   rm(hbeta, hy, lam, nX, objlasso, sigmanew);

   est$fitted.values = as.vector(X %*% est$coefficients);
   est$residuals = y - est$fitted.values;
   est$type = "regression";
   est$call = match.call();
   class(est) = "scalreg";
   #=================================================================
   est
}

bastah <- function(X, y, family = "gaussian", mcorr = "holm", N = 10000, ncores = 4, verbose = FALSE) {
   # Removing variables with a constant value by finding levels in the selected data
   n = nrow(X);
   p = ncol(X);
   numLevels = vector(mode = "numeric", p);
   for (i in 1:p) {
      numLevels[i] = length(unique(X[, i]));
   }
   selection = which(numLevels > 1);
   X = X[, selection];
   rm(numLevels);
   gc();

   n = nrow(X);
   p = ncol(X);

   if(verbose){
      print(paste(Sys.time(),"Pre-processing data"));
   }
   # Normalizing data using ridge regression
   if(family == "binomial"){
      X = scale(X, center = TRUE, scale = FALSE);
      D = rep(0, p);
      for(i in 1:p){
         D[i] = sqrt((n - 1) / (t(X[, i]) %*% X[, i]));
      }
      D = diag(D);
      X = X %*% D;
      rm(D);
      gc();
   } else if(family != "gaussian"){
      stop("Invalid family for response variable");
   }

   if(verbose){
      print(paste(Sys.time(),"Computing precision matrix"));
   }
   # Running BigQuic to calculate Z
   bqTheta.hat <- BigQuic(X, lambda = 0.99, use_ram = TRUE, numthreads = ncores);
   Theta.hat <- bqTheta.hat$precision_matrices[[1]];
   Z = Theta.hat %*% t(X);
   Z = as.matrix(t(Z));
   rm(bqTheta.hat, Theta.hat);
   gc();

   X <- apply(X, 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))

   if(family == "binomial"){
      fitnet = cv.glmnet(X, y, family = "binomial", standardize = FALSE);
      glmnetfit = fitnet$glmnet.fit;
      netlambda.min = fitnet$lambda.min;
      netpred = predict(glmnetfit, X, s = netlambda.min, type = "response");
      betahat = predict(glmnetfit, X, s = netlambda.min, type = "coefficients");
      betahat = as.vector(betahat);
      pihat = netpred[ , 1 ];
      diagW = pihat * (1 - pihat);
      W = diag(diagW);
      xl = cbind(rep(1, n), X);
      ## Adjusted design matrix
      X = sqrt(diagW) * X;
      ## Adjusted response
      y = sqrt(diagW) * (xl %*% betahat + solve( W, y - pihat));
      sigma = 1;
      rm(fitnet, betahat, diagW, glmnetfit, netlambda.min, netpred, pihat, W, xl);
      gc();
   } else {
      sigma = NULL;
   }

   ## center the columns the response to get rid of the intercept
   X <- scale(X, center = TRUE, scale = FALSE);
   y <- scale(y, scale = FALSE);
   y <- as.numeric(y);

   if(verbose){
      print(paste(Sys.time(),"Running lasso"));
   }
   # Bias estimate based on lasso estimate
   scaledlassofit = bastah.scalreg(X = X, y = y);
   betalasso <- scaledlassofit$coefficients

   ## Get estimated standard deviation
   if(is.null(sigma)){
      sigmahat = scaledlassofit$hsigma;
   } else {
      sigmahat = sigma;
   }
   rm(scaledlassofit);
   gc();

   # Rescale the Z appropriately such that such that t(Zj) Xj/n = 1 for all j
   scaleZ = vector(mode = "numeric", n);
   for(i in 1:ncol(X)){
      scaleZ[i] = (Z[, i] %*% X[, i]) / n;
   }
   Z = scale(Z, center = FALSE, scale = scaleZ);
   rm(scaleZ);
   gc();

   if(verbose){
      print(paste(Sys.time(),"Projection estimator and bias"));
   }
   # Projection estimator and bias
   bproj = t(Z) %*% y / n;

   if (ncores > 1 & !requireNamespace("doMC", quietly = TRUE)) {
      warning("doMC package is not available using single core");
      ncores = 1;
   }
   if(ncores > 1){
      doMC::registerDoMC(ncores);
      AA = crossprod(Z, X);
      bias = as.numeric( as.vector( foreach( j = 1:p ) %dopar%
      {
         ( as.vector( AA[ j, -j ] ) ) %*% betalasso[ -j ] / n
      } ) )
      rm(AA);
      gc();
   } else {
      bias = numeric(p);
      for(j in 1:p){
         bias[j] = (t(Z[ ,j]) %*% X[, -j]) %*% betalasso[-j] / n;
      }
   }

   bproj = bproj - bias;
   bproj = betalasso + bproj;
   rm(bias, betalasso);
   gc();

   # Determine normalization factor
   scaleb = n / (sigmahat * sqrt(colSums(Z ^ 2)));

   # Calculate p-value
   pval = 2 * pnorm(abs(bproj * scaleb), lower.tail = FALSE);
   rm(scaleb);

   # Multiple testing correction
   if(mcorr == "WY"){
      ## Westfall-Young like procedure as in ridge projection method,
      ## P.Buhlmann & L.Meier
      ## method related to the Westfall - Young procedure
      ## constants left out since we"ll rescale anyway
      ## otherwise cov2 <- crossprod(Z)/n
      pcorr <- p.adjust.wy( cov = crossprod(Z), pval = pval, N = N )
   } else if(mcorr %in% p.adjust.methods) {
      pcorr = p.adjust(pval, method = mcorr);
   } else {
      stop("Unknown multiple correction method specified");
   }
   rm(Z);
   gc();

   out = list(pval        = as.vector(pval),
              pval.corr   = pcorr,
              sigmahat    = sigmahat,
              bhat        = bproj,
              selection   = selection
   );
   names(out$pval) = names(out$pval.corr) = names(out$bhat) = names(out$betahat) = colnames(X);
   class(out) <- "bastah";
   if(verbose){
      print(paste(Sys.time(),"Done"));
   }

   return(out);
}

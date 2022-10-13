    # The code below is the implementation of three pricing methods:
    # (attachted paper provide detailed description)
    # a) Stotzel`s method (Stoetzel, 1970)
    # b) van Vestendorp`s method (Vestendorp, 1976)
    # c) Lipovetsky`s method (Lipovetsky, 2011)
    #
    # The data set is the survey data which was aimed on investigating
    # into price intervals of private tutor`s servives, n = 37
    # 
    # The code provides a final results of the use of these methods.
    #
    # Literature
    #
    # Stoetzel J. Psychological/ Sociological Aspects of Price. Bernard 
    # Taylor and Gordon Wills (Eds.).Pricing Research, Princeton, NJ: Brandon Systems 
    # Press, 1970. pp. 70ñ74
    #
    # Westendorp P.H. NSS ñ Price Sensitivity Meter (PSM) ñ A New 
    # Approach to Consumer Perception of Prices. Venice Congress Main Sessions, 
    # European Marketing Research Society (ESOMAR), Amsterdam, 1976. –. 139ñ167
    #
    # Lipovetsky, S., Magnan, S., & Zanetti-Polzi, A. (2011). Pricing models in 
    # marketing research. Intelligent Information Management. 
    # 
    # ‘¿Ã»À»ﬁ “”“ Õ¿œ»ÿ» — »Ã≈Õ≈Ã —¬Œ»Ã » ——€ÀÍÛ Ì‡ œÓ◊“” Œ¡€◊ÕŒ “¿  ƒ≈À¿ﬁ“
    
    library(ggplot2)        
    
    # Stoetzel method implementation
    
    Stoetzel_method <- function(X, Y)
    {
        N1 = length(X);
        N2 = length(Y);
        Ix <- rep(1, N1);
        Iy <- rep(2, N2);   
        Z <- c(X, Y);
        Iz <- c(Ix, Iy);    
        Z <- sort(Z, index.return = T);
        S <- data.frame(Z$x, Iz[Z$ix]);    
        names(S) <- c("Z", "Ind");
        dx <- 1 / N1;
        dy <- 1 / N2;
        Diff <- rep(0, N1 + N2);
        Diff[1] <- dx;
        #It is required to implement static or dynamic step for CDF due to precense of identical values
        
        for (i in seq(2, N1 + N2))
        {
            if (S$Ind[i] == 1)
            {
                Diff[i] <- Diff[i - 1] + dx;
            }
            else
            {
                Diff[i] <- Diff[i - 1] - dy;
            }
        }
        i = 1;
        N = length(Diff);
        while (i <= N)
        {
            j = i;
            while (S$Z[j] == S$Z[j + 1])
            {
                j = j + 1;
            }
            Diff[i] = Diff[j];
            if (S$Z[i] == S$Z[i + 1])
            {
                S$Z[seq(i + 1, N - (j - i))] <- S$Z[seq(j + 1, N)];
                Diff[seq(i + 1, N - (j - i))] <- Diff[seq(j + 1, N)];
                N = N - (j - i);
            }
            i = i + 1;
        }
        
        d1 <- NA;
        d2 <- NA;
        # Looking for MINIMAL solution d(p) = 0.5        
        for (i in seq(2, N))
        {
            if ((Diff[i - 1] < 0.5) && (Diff[i] > 0.5))
            {
                break;
            }
        }           
        if (i != N)
        {
            alpha = (Diff[i] - Diff[i - 1]) / (S$Z[i] - S$Z[i - 1]);
            beta = Diff[i] - alpha * S$Z[i];    
            d1 = (0.5 - beta) / alpha;        
        }        
        
        # Looking for MAXIMAL solution d(p) = 0.5
        for (i in seq(N, 1, by = -1))
        {
            if ((Diff[i - 1] > 0.5) && (Diff[i] < 0.5))
            {
                break;
            }
        }
        if (i != 1)
        {
            alpha = (Diff[i] - Diff[i - 1]) / (S$Z[i] - S$Z[i - 1]);
            beta = Diff[i] - alpha * S$Z[i];    
            d2 = (0.5 - beta) / alpha;        
        }            
        
        Stoetzel_method <- c(d1, d2)   
    }
    
    # Lipovetskiy`s method
    
    Lipovetskii_method <- function(data)
    {
        # Computing Pi, i = 1,4
        N <- nrow(data);
        dx = 1 / N;
        P1 = seq(dx, 1, by = dx); 
        P2 = seq(dx, 1, by = dx); 
        P3 = seq(dx, 1, by = dx); 
        P4 = seq(dx, 1, by = dx)
        data <- cbind(data, P1, P2, P3, P4)
        
        data$X1 <- sort(data$X1);
        data$X2 <- sort(data$X2);
        data$X3 <- sort(data$X3);
        data$X4 <- sort(data$X4);
        
        # Linear Modeling
        logP1 = log(data$P1 / (1 - data$P1)); logP2 = log(data$P2 / (1 - data$P2));
        logP3 = log(data$P3 / (1 - data$P3)); logP4 = log(data$P4 / (1 - data$P4));
        data <- cbind(data, logP1, logP2, logP3, logP4)
        data[data == Inf] <- NA
        data <- na.omit(data)  
        lm1 = lm(data = data, logP1 ~ log(X1))
        lm2 = lm(data = data, logP2 ~ log(X2))
        lm3 = lm(data = data, logP3 ~ log(X3))
        lm4 = lm(data = data, logP4 ~ log(X4))
        
        # Lipovetsky`s Algorithm
        
        a1_heat = lm1$coefficients[1]; b1_heat = lm1$coefficients[2];
        a2_heat = lm2$coefficients[1]; b2_heat = lm2$coefficients[2];
        a3_heat = lm3$coefficients[1]; b3_heat = lm3$coefficients[2];        
        a4_heat = lm4$coefficients[1]; b4_heat = lm4$coefficients[2];
        
        Xmin <- min(data$X1, data$X2, data$X3, data$X4)
        Xmax <- max(data$X1, data$X2, data$X3, data$X4)
        
        X <- seq(Xmin, Xmax, by = 0.1)  
        Y1 <- 1 / (1 + exp(-(a1_heat + b1_heat * log(X))));
        Y2 <- 1 / (1 + exp(-(a2_heat + b2_heat * log(X))));
        Y3 <- 1 / (1 + exp(-(a3_heat + b3_heat * log(X))));
        Y4 <- 1 / (1 + exp(-(a4_heat + b4_heat * log(X))));
        #d <- Y1 - Y4;
        d <- Y2 - Y3;
        
        P1 <- -1; P2 <- -1;
        for (i in seq(2, length(X)))
        {
            if (  (d[i - 1] < 0.5) && (d[i] > 0.5) )
                break;
        }
        
        b <- (d[i] - d[i - 1]) / (X[i] - X[i - 1]);
        a = d[i] - b * X[i];
        P1 <- (0.5 - a) / b;
        for (i in seq(length(X), 2, by = -1))
        {
            if (  (d[i - 1] > 0.5) && (d[i] < 0.5) )
                break;
        }
        b <- (d[i] - d[i - 1]) / (X[i] - X[i - 1]);
        a = d[i] - b * X[i];
        P2 <- (0.5 - a) / b;    
        Lipovetskii_method <- c(P1, P2);
    }
    
    # auxiliary function for PSM method
    
    Normalize <- function(X, Fx)
    {
        i = 1;
        N = length(X);
        dx <- 1/N;        
        while (i <= N)
        {
            j = i;
            while (X[j] == X[j + 1])
            {
                j = j + 1;
            }                        
            Fx[i] = Fx[j];
            if (X[i] == X[i + 1])
            {
                X[seq(i + 1, N - (j - i))] <- X[seq(j + 1, N)];
                Fx[seq(i + 1, N - (j - i))] <- Fx[seq(j + 1, N)];
                N = N - (j - i);
            }
            i = i + 1;
        }        
        Normalize <- data.frame(X[seq(1, N)], Fx[seq(1, N)]);
    }
    
    # Another auxiliary function for PSM method
    
    Intersect <- function(X, Fx, Y, Fy)    
    {        
        i <- 2; j <- 2;
        N <- length(X); 
        M <- length(Y);        
        
        for (i in seq(2, N))
        {
            for (j in seq(2, M))
            {
                a1 <- (Fx[i] - Fx[i - 1]) / (X[i] - X[i - 1]);
                b1 <- Fx[i] - a1 * X[i];
                
                a2 <- (Fy[j] - Fy[j - 1]) / (Y[j] - Y[j - 1]);
                b2 <- Fy[j] - a2 * Y[j];
                
                x = (b2 - b1) / (a1 - a2);
                y = a1 * x + b1;        
                
                if ((x <= X[i]) && (x >= X[i - 1]) && (x <= Y[j]) && (x >= Y[j - 1]))
                {
                    X_inter <- x;
                    Y_inter <- y;
                }
            }
        }
        
        Intersect <- c(X_inter, Y_inter);
    }
    
    PSM_method <- function(data)
    {
        X1 <- data[[1]];
        X2 <- data[[2]];
        X3 <- data[[3]];
        X4 <- data[[4]];
        X1 <- sort(X1); X2 <- sort(X2);
        X3 <- sort(X3); X4 <- sort(X4);
        
        N <- length(X1);    
        Fx <- seq(1/N, 1, 1/N);                
        
        D <- Normalize(X1, Fx);
        X1 <- D$X; S1 <- 1 - D$Fx;    
        D <- Normalize(X2, Fx);
        X2 <- D$X; S2 <- 1 - D$Fx;    
        D <- Normalize(X3, Fx);    
        X3 <- D$X; F3 <- D$Fx;    
        D <- Normalize(X4, Fx);    
        X4 <- D$X; F4 <- D$Fx;    
        
        S1 <- c(1, S1); X1 <- c(0, X1);
        S2 <- c(1, S2); X2 <- c(0, X2);
        F3 <- c(0, F3); X3 <- c(0, X3);
        F4 <- c(0, F4); X4 <- c(0, X4);
        
                
        # Point of minimal cost
        PMC <- Intersect(X1, S1, X3, F3)[1];
        # Point of maximal espensiveness
        PME <- Intersect(X2, S2, X4, F4)[1];
        # Indifference price point
        IPP <- Intersect(X2, S2, X3, F3)[1];
        # Optimal price point
        OPP <- Intersect(X1, S1, X4, F4)[1];
        
        PSM_method = c(PMC, PME, IPP, OPP)
    }
    
    ############################################################################
    ######################## TESTING THE METHODS ON THE DATA SET ###############
    ############################################################################
    
    #setwd(":\\")
    data <- read.csv("tutor_price_data.csv", head = T, sep = ";")    
    names(data) <- c("X1", "X2", "X3", "X4")        
    head(data)
    
    ST_res <- Stoetzel_method(data[[2]], data[[3]])    
    ST_res
        
    LIP_res <- Lipovetskii_method(data)
    LIP_res
    
    PSM_res <- PSM_method(data)
    PSM_res
    

Ln <- function(Tn=Tn,n=n)
{
    B <- dim(Tn)[1]
    r <- expand.grid(rep(list(c(0,1)), dim(Tn)[2])) # 2^p possible a
    t <- apply(r,1,sum)
    r <- r[-which(t==0|t==dim(Tn)[2]),]
    Ln <- matrix(NA,ncol=1,nrow=B)


    for (b in 1:B)
        {
            gg <- c()
            
            Sna_aap1 <- cov(Tn)
            
            for (i in 1:dim(r)[1])
                {
                    tmp <- which(r[i,]==1)
                    if (length(tmp)>0 & length(tmp)<dim(Tn)[2])
                        {
                            Tna <- matrix(Tn[,tmp],ncol=length(tmp))
                            Tnap <- matrix(Tn[,-tmp],ncol=(dim(Tn)[2]-length(tmp)))
                            Snaa <- cov(Tna) # Snaa
                            Snaap <- cov(Tna,Tnap) # Snaa'
                            Snapa <- t(Snaap)  # Sna'a
                            Snapap <- cov(Tnap) # Sna'a'
                            
                            Sna_aap <- Snaa-Snaap %*% solve(Snapap) %*% Snapa
                        
                        }
                    
                    Tna_ap <- Tna[b,] - Snaap %*% solve(Snapap) %*% Tnap[b,]
                    
                    tt <- ifelse(Tna_ap>0,1,0)
                    vv <- ifelse(solve(Snapap) %*% Tnap[b,]<=0,1,0)
                    Ina <- ifelse((sum(tt)+sum(vv))==(length(tt)+length(vv)),1,0)
                    
                    gg <- append(gg,(n*n)*Ina * t(Tna_ap) %*% solve(Sna_aap) %*% Tna_ap) # incluir n^2
                }
            
            # case where all a's = 1
            Tna_ap1 <- Tn[b,]
            tt1 <- ifelse(Tna_ap1>0,1,0)
            Ina1 <- ifelse(sum(tt1)==length(tt1),1,0)
            
            Ln[b] <- sum(gg) + (n*n)*Ina1 * t(Tna_ap1) %*% solve(Sna_aap1) %*% Tna_ap1
        }
    return(Ln)
}

    
    
    

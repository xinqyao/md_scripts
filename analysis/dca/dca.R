#' Direct Coupling Analysis (DCA)
#'
#' Perform direct-coupling analysis on a multiple sequence alignment
#'
#' This function is converted from the original Matlab file, dca.m, with slight
#' modificiations. Following are comments from the original file:
#'
#' Direct Coupling Analysis (DCA)
#'
#' function dca(inputfile , outputfile)
#' 
#' INPUTS: 
#'   inputfile  - file containing the FASTA alignment
#'   outputfile - file for dca results. The file is composed by N(N-1)/2 
#'                (N = length of the sequences) rows and 4 columns: 
#'                residue i (column 1), residue j (column 2),
#'                MI(i,j) (Mutual Information between i and j), and 
#'                DI(i,j) (Direct Information between i and j).
#'                Note: all insert columns are removed from the alignment.
#'
#' SOME RELEVANT VARIABLES:
#'   N        number of residues in each sequence (no insert)
#'   M        number of sequences in the alignment
#'   Meff     effective number of sequences after reweighting
#'   q        equal to 21 (20 aminoacids + 1 gap)
#'   align    M x N matrix containing the alignmnent
#'   Pij_true N x N x q x q matrix containing the reweigthed frequency
#'            counts.
#'   Pij      N x N x q x q matrix containing the reweighted frequency 
#'            counts with pseudo counts.
#'   C        N(q-1) x N(q-1) matrix containing the covariance matrix.
#'
#'
#' Copyright for this implementation: 
#'             2011/12 - Andrea Pagnani and Martin Weigt
#'                       andrea.pagnani@gmail.com 
#'                       martin.weigt@upmc.fr
#' 
#' Permission is granted for anyone to copy, use, or modify this
#' software and accompanying documents for any uncommercial
#' purposes, provided this copyright notice is retained, and note is
#' made of any changes that have been made. This software and
#' documents are distributed without any warranty, express or
#' implied. All use is entirely at the user's own risk.
#'
#' Any publication resulting from applications of DCA should cite:
#'
#'     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander, 
#'     R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
#'     analysis of residue co-evolution captures native contacts across 
#'     many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
#'
#'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'
#' @param inputfile file containing the FASTA alignment (e.g., from the "full"
#'        alignement of a PFAM famility; Note that an alignment file having 
#'        mixed "." and "-" for gaps, as well small cases for insertions, is required).
#' @param outputfile file for dca results.
#' @param outputfile2 file for Pi results (Optional).
#'
#' @return
#'
#' @references
#'        Morcos, F. et al. (2011) \emph{PNAS} \bold{108}, E1293-1301.
#'
#' @examples
#'
#' @keywords analysis 
dca <- function (inputfile, outputfile, outputfile2=NULL) {

    require(bio3d)

    pseudocount_weight = 0.5  # relative weight of pseudo count   
    theta = 0.2               # threshold for sequence id in reweighting

    
    out <- .return_alignment(inputfile)
    N <- out$N; M <- out$M; q <- out$q; align <- out$align
    out <- .Compute_True_Frequencies(align,M,N,q, theta)
    Pij_true <- out$Pij_true; Pi_true <- out$Pi_true; Meff <- out$Meff
    
    cat(sprintf('### N = %d M = %d Meff = %.2f q = %d\n', N,M,Meff,q))
    out <- .with_pc(Pij_true,Pi_true,pseudocount_weight,N,q)
    Pij <- out$Pij; Pi <- out$Pi
    C <- .Compute_C(Pij,Pi,N,q)
#    invC <- inv(C)
    invC <- solve(C)
    
    di <- .Compute_Results(Pij, Pi,Pij_true, Pi_true, invC, N, q)
#    write.table(di, file=outputfile, row.names=FALSE, col.names=FALSE)
    cat("", file=outputfile)
    for(i in 1:nrow(di)) {
       cat(sprintf("%d %d %g %g\n", di[i, 1], di[i, 2], di[i, 3], di[i, 4]), file=outputfile, append=TRUE)
    }
    if(!is.null(outputfile2)) {
#       write.table(Pi, file=outputfile2, row.names=FALSE, col.names=FALSE)
       cat("", file=outputfile2)
       for(i in 1:nrow(Pi)) { 
          for(j in 1:ncol(Pi)) {
             cat(sprintf("%g ", Pi[i, j]), file=outputfile2, append=TRUE)
          }
          cat("\n", file=outputfile2, append=TRUE)
       }
    }
    return( list(di=di, Pi=Pi) )
}

.return_alignment <- function(inputfile) {
# reads alignment from inputfile, removes inserts and converts into numbers

    align_full <- bio3d::read.fasta(inputfile, rm.dup=FALSE, to.upper=FALSE, to.dash=FALSE)
    M <- length(align_full$id)
    ind <- align_full$ali[1, ] %in% c("-", LETTERS)
    N <- sum(ind)
    Z <- matrix(0, nrow=M, ncol=N)

    for(i in 1:M) {
        counter = 0
        for(j in 1:length(ind)) {
            if( ind[j] ) {
                counter <- counter + 1
                Z[i, counter] <- .letter2number( align_full$ali[i, j] )
            }
        }
    }
    q <- max(Z)

    return(list(N=N, M=M, q=q, align=Z))
}

.Compute_Results <- function(Pij,Pi,Pij_true,Pi_true,invC, N,q) {
# computes and prints the mutual and direct informations

    di <- NULL
    for (i in 1:(N-1)) {
        for (j in (i+1):N) {
            # mutual information
            out <- .calculate_mi(i,j,Pij_true,Pi_true,q)
            MI_true <- out$MI_true; si_true <- out$si_true 
            sj_true <- out$sj_true
            
            # direct information from mean-field
            W_mf <- .ReturnW(invC,i,j,q)
            DI_mf_pc <- .bp_link(i,j,W_mf,Pi,q)
            di <- rbind(di, data.frame(i=i, j=j, MI=MI_true, DI=DI_mf_pc))
        }
    }
    return( di )
}

.Compute_True_Frequencies <- function(align,M,N,q,theta) {
# computes reweighted frequency counts

    W <- rep(1,M)
    if( theta > 0.0 ) {
#        W <- (1./(1+colSums(apply(align, 1, function(x) colSums(x != t(align)) / N)<theta)))
        W <- (1./(colSums(apply(align, 1, function(x) colSums(x != t(align)) / N)<theta)))
    }
    Meff <- sum(W)

    Pij_true <- array(0, dim=c(N,N,q,q))
    Pi_true <- matrix(0, N, q)

    for (j in 1:M) {
        for (i in 1:N) {
            Pi_true[i,align[j,i]] <- Pi_true[i,align[j,i]] + W[j]
        }
    }
    Pi_true <- Pi_true/Meff

    for (l in 1:M) {
        for (i in 1:(N-1)) {
            for (j in (i+1):N) {
                Pij_true[i,j,align[l,i],align[l,j]] <- Pij_true[i,j,align[l,i],align[l,j]] + W[l]
                Pij_true[j,i,align[l,j],align[l,i]] <- Pij_true[i,j,align[l,i],align[l,j]]
            }
        }
    }
    Pij_true <- Pij_true/Meff

    scra <- diag(q)
    for (i in 1:N) {
        for (alpha in 1:q) {
            for (beta in 1:q) {
                Pij_true[i,i,alpha,beta] <- Pi_true[i,alpha] * scra[alpha,beta]
            }
        }
    }
    return( list(Pij_true=Pij_true, Pi_true=Pi_true, Meff=Meff) )
}

.letter2number <- function(a) {
    # full AA alphabet
    switch(a,
        '-'=1,
        'A'=2,    
        'C'=3,
        'D'=4,
        'E'=5,
        'F'=6,
        'G'=7,
        'H'=8,
        'I'=9,
        'K'=10,
        'L'=11,
        'M'=12,
        'N'=13,
        'P'=14,
        'Q'=15,
        'R'=16,
        'S'=17,
        'T'=18,
        'V'=19,
        'W'=20,
        'Y'=21,
        1)
}

.with_pc <- function(Pij_true, Pi_true, pseudocount_weight,N,q) {
# adds pseudocount

    Pij <- (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*array(1, dim=c(N,N,q,q))
    Pi <- (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*matrix(1,N,q)

    scra <- diag(q)

    for (i in 1:N) {
        for (alpha in 1:q) {
            for (beta in 1:q) {
               Pij[i,i,alpha,beta] <-  (1.-pseudocount_weight)*Pij_true[i,i,alpha,beta] + pseudocount_weight/q*scra[alpha,beta]
            }
        }
    }
    return( list(Pij=Pij, Pi=Pi) )
}

.Compute_C <- function(Pij,Pi,N,q) {
# computes correlation matrix

    C <- matrix(0, nrow=N*(q-1), ncol=N*(q-1))
    for (i in 1:N) {
        for (j in 1:N) {
            for (alpha in 1:(q-1)) {
                for (beta in 1:(q-1)) {
                     C[.mapkey(i,alpha,q),.mapkey(j,beta,q)] <- Pij[i,j,alpha,beta] - Pi[i,alpha]*Pi[j,beta]
                }
            }
        }
    }
    return( C )
}

.mapkey <- function(i,alpha,q) {
    (q-1)*(i-1)+alpha
}

.calculate_mi <- function(i,j,P2,P1,q) {
# computes mutual information between columns i and j

    M <- 0.
    for (alpha in 1:q) {
        for (beta in 1:q) {
             if( P2[i,j,alpha,beta]>0 ) {
                M <- M + P2[i,j,alpha, beta]*log(P2[i,j, alpha, beta] / P1[i,alpha]/P1[j,beta])
            }
        }
    }

    s1 <- 0.
    s2 <- 0.
    for (alpha in 1:q) {
        if( P1[i,alpha]>0 ) {
            s1 <- s1 - P1[i,alpha] * log(P1[i,alpha])
        }
        if( P1[j,alpha]>0 ) {
            s2 <- s2 - P1[j,alpha] * log(P1[j,alpha])
        }
    }
    return( list(MI_true=M, si_true=s1, sj_true=s2) )
}

.ReturnW <- function(C,i,j,q) {
# extracts coupling matrix for columns i and j

    W <- matrix(1,q,q)
    W[1:(q-1),1:(q-1)] <- exp( -C[.mapkey(i,1:(q-1),q),.mapkey(j,1:(q-1),q)] )
    W
}

.bp_link <- function(i,j,W,P1,q) {
# computes direct information

    out <- .compute_mu(i,j,W,P1,q)
    mu1 <- out$mu1; mu2 <- out$mu2
    DI <- .compute_di(i,j,W, mu1,mu2,P1)
    return (DI)
}

.compute_mu <- function(i,j,W,P1,q) {

    epsilon=1e-4
    diff =1.0
    mu1 = matrix(1,1,q)/q
    mu2 = matrix(1,1,q)/q
    pi = P1[i,, drop=FALSE]
    pj = P1[j,, drop=FALSE]

    while ( diff > epsilon ) {

        scra1 = mu2 %*% t(W)
        scra2 = mu1 %*% W

        new1 = pi/scra1
        new1 = new1/sum(new1)

        new2 = pj/scra2
        new2 = new2/sum(new2)

        diff = max( max( abs( new1-mu1 ), abs( new2-mu2 ) ) )

        mu1 = new1
        mu2 = new2

    }
    return( list(mu1=mu1, mu2=mu2) )
}

.compute_di <- function(i,j,W, mu1,mu2, Pia) {
# computes direct information

    tiny = 1.0e-100

    Pdir = W * (t(mu1) %*% mu2)
    Pdir = Pdir / sum(Pdir)

    Pfac = t(Pia[i,,drop=FALSE]) %*% Pia[j,,drop=FALSE]

    DI = sum(diag( t(Pdir) %*% log( (Pdir+tiny)/(Pfac+tiny) ) ))

}

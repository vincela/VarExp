###############################################################################
######## Estimate the genotype correlation matrix from reference panel ########
###############################################################################

#' Code additively the genotypes
#' 
#' Sum the two columns per individual of the input matrix (each corresponding 
#' to one allele) to get the additive coding of the individual at the locus.
#'
#' @param genoMat is the matrix of genotype with two colums per individual
#' Rows are variants and columns are individuals
#' @param lind is the list of individual in the genotype matrix
#'
#' @return The additively coded genotype matrix
#'
#' @examples
#' genoMat = matrix(rbinom(5*10,1,runif(10,0,1)),nrow=5,ncol=10)
#' lind = c("ind1","ind2","ind3","ind4","ind5")
#' colnames(genoMat) = rep(lind,each=2)
#' addCodedMat = getAdditivelyCodedMatrix(genoMat, lind)
#'
getAdditivelyCodedMatrix = function(genoMat, lind) {
  res = matrix(0,nrow = dim(genoMat)[1],ncol=length(lind))
  getGeno = function(i,genoMat,lind,res) {
    res[,i] = apply(genoMat[,grepl(lind[i], colnames(genoMat))],1,sum)
  }
  res = mapply(getGeno,1:length(lind),MoreArgs=list(genoMat=genoMat,lind=lind,res=res))
  return(res)
}

#' Align the coding with the reference panel
#' 
#' Change the coding of the variant if the reference allele differs between 
#' the data and the reference panel. If reference alleles differs, 
#' switch 0s and 2s. Otherwise, no change.
#'
#' @param x is the matrix of additively coded genotypes
#' Rows are variants and columns are individuals
#' @param v is a vector of boolean indicating wether reference alleles differ or not
#'
#' @return The additively coded genotyped matrix with coded allele in the matrix
#' corresponding to the reference allele in the reference panel
#'
#' @examples
#' genoMat = matrix(rbinom(5*5,2,runif(5,0,1)),nrow=5,ncol=5)
#' sameRefAllele = rbinom(5,1,0.7)
#' newGenoMat = changeCoding(genoMat, sameRefAllele)
#'
changeCoding = function(x, v) {
  tt = x
  for (i in 1:dim(x)[1]) {
    if (!v[i]) {
      tt[i,] = -tt[i,] + 2
    }
  }
  return(tt)
}

#' Compute the genotype correlation matrix.
#'
#' From a set of variants identified by the pair (chromosome, position), extract 
#' the SNPs from a reference panel in the specified population.
#'
#' Currently, this is hard-coded to access 1000 Genomes phase3 data hosted by
#' Brian Browning (author of BEAGLE):
#'
#' \url{http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/}
#'
#' This implementation discards multi-allelic markers that have a "," in the
#' ALT column.
#'
#' Position must be given in GRCh37 genome build.
#'
#' The \code{pop} can be any of: ACB, ASW, BEB, CDX, CEU, CHB, CHS, CLM, ESN,
#' FIN, GBR, GIH, GWD, IBS, ITU, JPT, KHV, LWK, MSL, MXL, PEL, PJL, PUR, STU,
#' TSI, YRI. It can also be any super-population: AFR, AMR, EAS, EUR, SAS.
#'
#' Then, code additively the genotype and modify the additively coded allele if 
#' reference alleles differ between data and reference panel
#' and finally compute the correlation matrix.
#'
#' Physical position must in
#'
#' @param lchr is a vector with the chromosome number of the variants to extract
#' @param lpos is a vector with the physical position of the variants to extract
#' @param lrefall is a vector with the reference allele of the variant in the data.
#' @param pop is the 1000 Genomes code of the population in which data must be extracted
#'
#' @return The genotype correlation matrix of the specified variants
#'
#' @examples
#' chrom = c(8,4)
#' phys_pos = c(11843758,951947)
#' refall = c("A","T")
#' cor_matrix = getGenoCorMatrix(lchr = chrom,lpos = phys_pos , lrefall = refall, pop = "EUR")
#'
#' @export
#' 
getGenoCorMatrix = function(lchr,lpos, lrefall, pop) {
  if (length(lchr)>1) {
    lpop = rep(pop,length(lchr))
    referencedata = mapply(get_vcf,lchr,lpos,lpos,lpop)
    genoMat = Reduce("rbind", lapply(1:dim(referencedata)[2], 
                                     function(i) referencedata[,i]$geno))
    genoMap = Reduce("rbind", lapply(1:dim(referencedata)[2], 
                                     function(i) referencedata[,i]$meta))

    # Set allele to "t" when read allele is "TRUE"
    genoMap$REF = vapply(genoMap$REF, function(x) 
      if (!(as.character(x) %in% c("A","C","G"))) {return("t")} else {return(x)},"")

    #Keep only founders
    lind = referencedata[,1]$ind
    lind = referencedata[,1]$ind[which(lind$Paternal.ID == 0 & lind$Maternal.ID == 0),2]

    testMat = getAdditivelyCodedMatrix(genoMat,lind)
    testMat = changeCoding(testMat,tolower(as.character(lrefall)) == tolower(genoMap$REF))
    colnames(testMat) = lind
    return(stats::cor(t(testMat)))
  }
  else {
    return(1)
  }
}

#' Perform singular Value Decomposition on the correlation matrix
#'
#' @param cormat is the correlation matrix
#' @param k is the number of eigenvectors to keep. 
#'   Default is the correlation matrix rank.
#'
#' @return A list with :
#' \describe{
#'   \item{eigval}{A vector of the top \code{k} eigenvalues}
#'   \item{eigvev}{A matrix of the top \code{k} eigenvectors}
#' }
#'
#'
#' @examples
#' a = matrix(runif(10*10,1,5),nrow=10)
#' matcor = cor(a)
#' getMatCorSVD(matcor)
#' getMatCorSVD(matcor, k = 4)
#' 
getMatCorSVD = function(cormat, k = qr(cormat)$rank) {
  cormat.svd = svd(cormat, nu = 0, nv = k)
  list(eigval = cormat.svd$d[1:k], eigvec = cormat.svd$v)
}

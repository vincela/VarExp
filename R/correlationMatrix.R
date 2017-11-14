###################################################################################
############ Estimate genotype correlation matrix from reference panel ############
###################################################################################

#' Code additively the genotypes
#' Sum the two columns per individual of the input matrix (each corresponding to one chromosome)
#' to get the additive coding of the individual at the locus.
#'
#' @param genoMat is the matrix of genotype with two colums per individual
#' Rows are variants and columns are individuals
#' @param lind is the list of individual in the genotype matrix
#'
#' @return The additively coded genotyped matrix
#'
#' @examples
#' genoMat = matrix(rbinom(5*10,1,runif(10,0,1)),nrow=5,ncol=10)
#' lind = c("ind1","ind2","ind3","ind4","ind5")
#' colnames(genoMat) = rep(lind,each=2)
#' addCodedMat = getAdditivelyCodedMatrix(genoMat, lind)
#'
#'@export

getAdditivelyCodedMatrix = function(genoMat,lind) {
  res = matrix(0,nrow = dim(genoMat)[1],ncol=length(lind))
  getGeno = function(i,genoMat,lind,res) {
    res[,i] = apply(genoMat[,grepl(lind[i], colnames(genoMat))],1,sum)
  }
  res = mapply(getGeno,1:length(lind),MoreArgs=list(genoMat=genoMat,lind=lind,res=res))
  return(res)
}

#' Change the coding of the variant
#' if the reference allele differs between the data and the reference panel
#'
#' If reference alleles differs, new coding = 2 - old_coding
#' Otherwise, no change
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
#'@export

changeCoding = function(x,v) {
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
#' From a set of variants identified by the pair (chromosome, position), extract the SNPs from a reference
#' panel in the specified population.
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
#' Then, code additively the genotype and modify the additively coded allele if reference alleles
#' differ between data and reference panel
#' and finally compute the correlation matrix.
#'
#' Physical position must in
#'
#' @param lrsid is a vector with the rs identifiers of the variants to extract
#' @param lchr is a vector with the chromosome number of the variants to extract
#' @param lpos is a vector with the physical position of the variants to extract
#' @param lrefall is a vector with the reference allele of the variant in the data.
#' @param pop is the 1000 Genomes code of the population in which data must be extracted
#' @param path path to the reference genotype file if locally downloaded
#' @param web logical indicating whether to use web access or not
#'
#' @return The genotype correlation matrix of the specified variants
#'
#' @examples
#'
#' ids = c("rs4841662", "rs34311866")
#' chrom = c(8,4)
#' phys_pos = c(11843758,951947)
#' refall = c("A","T")
#' cor_matrix = getGenoCorMatrix(lrsid = ids, lchr = chrom,lpos = phys_pos , lrefall = refall, pop = "EUR")
#'
#' @export

getGenoCorMatrix = function(lrsid, lchr, lpos, lrefall, pop, path = "", web = TRUE) {
  if (length(lchr)>1) {
    lpop = rep(pop,length(lchr))
    referencedata = mapply(get_vcf, lchr, lpos, lpos, lpop, path, web)
    genoMat = Reduce("rbind", lapply(1:dim(referencedata)[2], function(i) referencedata[,i]$geno))
    genoMap = Reduce("rbind",lapply(1:dim(referencedata)[2], function(i) referencedata[,i]$meta))

    # Set allele to "t" when read allele is "TRUE"
    genoMap$REF = vapply(genoMap$REF, function(x) if (!(as.character(x) %in% c("A","C","G"))) {return("t")} else {return(x)},"")

    #Keep only founders
    lind = referencedata[,1]$ind
    lind = referencedata[,1]$ind[which(lind$Paternal.ID == 0 & lind$Maternal.ID == 0),2]

    testMat = getAdditivelyCodedMatrix(genoMat,lind)
    to_remove <- NULL
    for (i in 1:nrow(testMat)) {
      if (!(rownames(testMat)[i] %in% lrsid)) {
        to_remove   <- c(to_remove, i)
      }
    }
    if (!is.null(to_remove)) {
      testMat   <- testMat[-to_remove,]
      genoMap   <- genoMap[-to_remove,]
    }
    testMat = changeCoding(testMat,tolower(as.character(lrefall)) == tolower(genoMap$REF))
    colnames(testMat) = lind
    return(cor(t(testMat)))
  }
  else {
    return(1)
  }
}

#' Perform singular Value Decomposition on the correlation matrix
#'
#' @param cormat is the correlation matrix
#' @param k is the number of eigen vectors to keep. Default is the correlation matrix rank.
#'
#' @return A list with :
#' \describe{
#'   \item{eigval}{A vector of the top \code{k} eigenvalues}
#'   \item{eigvev}{A matrix of the top \code{k} eigen vectors}
#' }
#'
#'
#' @examples
#' a = matrix(runif(10*10,1,5),nrow=10)
#' matcor = cor(a)
#' getMatCorSVD(matcor)
#' getMatCorSVD(matcor, k = 4)
#'
#' @export

getMatCorSVD = function(cormat, k = qr(cormat)$rank) {
  cormat.svd = svd(cormat)
  l = list()
  l$eigval = cormat.svd$d[1:k]
  l$eigvec = cormat.svd$v[,1:k]
  return(l)
}

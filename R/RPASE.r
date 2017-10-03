#' Computing binomial normal marginal density using 32-points gaussian quadrature.
#' @source Pawitan Y. In all likelihood: statistical modelling and inference using likelihood[M]. Oxford University Press, 2001.
#' @param x a quantile
#' @param n number of trials
#' @param p probability of success on each trial
#' @param sig overdispersion
#' @return The marginal probability density of binomial normal distribution

dbinorm32 <- function(x,n,p,sig){
  xi<- c(.0483076656,.1444719615,.2392873622,.3318686022,
         .4213512761,.5068999089,.5877157572,.6630442669,
         .7321821187,.7944837959,.8493676137,.8963211557,
         .9349060759,.9647622555,.9856115115,.9972638618)
  xi<- c(-xi[16:1],xi)
  wi<- c(.0965400885,.0956387200,.0938443990,.0911738786,
          .0876520930,.0833119242,.0781938957,.0723457941,
          .0658222227,.0586840934,.0509980592,.0428358980,
          .0342738629,.0253920653,.0162743947,.0070186100)
  wi<- c(wi[16:1],wi)
  th <- log(p/(1-p))
  u<- xi*4*sig
  pu <- exp(th+u)/(1+exp(th+u))
  pu<- ifelse(pu<exp(-25),exp(-25),pu)
  pu<- ifelse(pu>1-exp(-25),1-exp(-25),pu)
  probu<- stats::dbinom(x,n,pu)*stats::dnorm(u,0,sig)
  4*sig*sum(wi*probu)
}

#' Create an unfiltered RPASE phased block list from VCF data and an annotation table
#' @param phased_vcf_data A data.frame with columns ordered according to the VCF standard containing information about phased SNPs
#' @param annotation A data.frame containing gene annotations ordered as a four-column BED format
#' @return A list of data.frames, each describing a phased block as determined from inputs
#' @details For further details, see the documentation for \code{\link{createPhasedBlockList}}

createUnfilteredPhasedBlockList <- function(phased_vcf_data,annotation){
  wm.1 = phased_vcf_data[,c(1,2,10)]
  wm.2= wm.1[wm.1[,3]!="./.",]

  a.2=gsub(".*1\\|0.*", "-1", wm.2[,3])
  a.3=gsub(".*0\\|1.*", "1", a.2)
  a.4=gsub(".*0/0.*", "2", a.3)
  a.5=gsub(".*1/1.*", "2", a.4)
  a.6=gsub(".*2/2.*", "2", a.5)
  a.7=gsub(".*3/3.*", "2", a.6)
  a.8=gsub(".*/.*", "0", a.7)
  a.9=gsub(".*\\|.*", "2", a.8)
  wm.data=data.frame(scaffold=wm.2[,1], position=wm.2[,2], index=as.numeric(a.9))
  wm.data=wm.data[wm.data$index!=2,]

  a=list()
  for (i in 1:nrow(annotation))
  {
    scaffold=as.character(annotation[i,1])
    vcf.sca=wm.data[wm.data$scaffold==scaffold,]
    star=annotation[i,2]
    end=annotation[i,3]
    a[[i]]=vcf.sca[((star<vcf.sca$position)&(vcf.sca$position<end)),]
  }
  names(a)=as.character(annotation[,4])

  y=list()
  for (i in names(a))
  {
    if(nrow(a[[i]])!=0)
    {
      y[[i]]=a[[i]]
    }
  }
  return(y)
}

#' Filter an RPASE phased block list
#' @param unfiltered_phased_block_list A list of phased blocks from createUnfilteredPhasedBlockList
#' @param min_phased_sites minimum number of phased SNPs to include within a phased block
#' @param min_coverage minimum coverage for each phased SNP included within a phased block
#' @param heter minimum minor-allele coverage for a phased SNP to be included within a phased block
#' @param phased_vcf_data A data.frame with columns ordered according to the VCF standard containing information about phased SNPs
#' @return A filtered RPASE phased block list
#' @details For further details, see  \code{\link{createPhasedBlockList}}

filterPhasedBlockList <- function(unfiltered_phased_block_list, min_phased_sites, min_coverage, heter, phased_vcf_data){
  new=list()
  for (q in 1:length(unfiltered_phased_block_list))
  {
    try.split.accoring.index=split(unfiltered_phased_block_list[[q]], cumsum(unfiltered_phased_block_list[[q]][,3] == 0))
    a=list()
    j=1
    for (i in 1: length(try.split.accoring.index))
    {
      if(sum(abs(try.split.accoring.index[[i]][try.split.accoring.index[[i]]$index!=2,]$index))>=min_phased_sites)
      {
        a[[j]]=try.split.accoring.index[[i]]
        j=j+1
      }
    }
    new[[q]]=a
  }
  names(new)=names(unfiltered_phased_block_list)

  new.without.null=list()
  for (p in names(new))
  {
    if(length(new[[p]])!=0)
    {
      new.without.null[[p]]=new[[p]]
    }
  }

  AD=function(u)
  {
    x.1=unlist(strsplit(u,':'))
    new=x.1[grep("^.*,.*$",x.1)]
    if(length(new)!=3)
    {
      new.3="NA"
    }else{
      new.2=as.numeric(unlist(strsplit(new[1],",",fixed=T)))
      if(length(new.2)!=3)
      {
        new.3="NA"
      }else{
        if(new.2[3]!=0)
        {
          new.3="NA"
        }else{
          new.3=paste(new.2[1],new.2[2],sep=",")
        }
      }
    }
    return (new.3)
  }
  wm.data=data.frame(scaffold=phased_vcf_data[,1], position=phased_vcf_data[,2], AD=sapply(as.character(phased_vcf_data[,10]), AD))
  wm.data.without.NA=wm.data[wm.data$AD!="NA",]
  aa=data.frame(do.call('rbind', strsplit(as.character(wm.data.without.NA$AD), ',', fixed=TRUE)))
  wm.data.AD=data.frame(scaffold=wm.data.without.NA$scaffold, position=wm.data.without.NA$position,
                        AD1=as.numeric(as.character(aa[,1])),AD2=as.numeric(as.character(aa[,2])))
  wm.data.AD.each.2=wm.data.AD[(wm.data.AD$AD1>heter)&(wm.data.AD$AD2>heter),]
  wm.data.AD.10=wm.data.AD.each.2[(wm.data.AD.each.2$AD1+wm.data.AD.each.2$AD2)>=min_coverage,]
  name.list=function(w)
  {
    length=length(w)
    names(w)=c(1:length)
    return(w)
  }
  b=sapply(new.without.null, name.list)
  if(all(sapply(b,length)==3)){
    b.2=do.call(rbind,b)
    b.3=data.frame(b.2, gene=rep(names(b),sapply(b,nrow)))
    rownames(b.3)=NULL
  }else{
    b.1=unlist(b,recursive=FALSE)
    b.2=do.call(rbind,b.1)
    b.3=data.frame(b.2, gene=rep(names(b.1),sapply(b.1,nrow)))
    rownames(b.3)=NULL
  }
  c=stats::na.omit(merge(b.3, wm.data.AD.10, by=c("scaffold","position"),all=TRUE))
  d=c[c$index==1|c$index==-1,]

  duo.1=split(d,d$gene,drop=T)

  twist.phase=function(x)
  {
    if(length(unique(x$index)) == 1)
    {
      return (x)
    }else{
      twist=x[which(x$index==-1),]
      twist.1=twist[,c(5,6)]
      twist[,5]=twist.1[,2]
      twist[,6]=twist.1[,1]
      twist.2=rbind(twist,x[which(x$index!=-1),])
      twist.3=twist.2[order(twist.2$position),]
      return (twist.3)
    }
  }
  twist.counts.list=lapply(duo.1,twist.phase)
  return(twist.counts.list)
}

#' Create a filtered RPASE phased block list from VCF data and an annotation table
#' @param phased_vcf_data A data.frame with columns ordered according to the VCF standard containing information about phased SNPs for a single sample
#' @param annotation A data.frame containing annotations with each line describing a study unit, typically a gene, ordered as a four-column BED format
#' @param min_phased_sites minimum number of phased SNPs to include within a phased block
#' @param min_coverage minimum coverage for each phased SNP included within a phased block
#' @param heter minimum minor-allele coverage for a phased SNP to be included within a phased block as heterzugous SNP
#' @return A list of data.frames, each describing a phased block, typically a gene. Each dataframe contains 6 columns: chromosome, position, index of phasing, gene, AD1 (allele depth of first allele), and AD2 (allele depth of second allele).
#' @details phased_vcf_data describes a single sample, and should thus have 10 columns corresponding to VCF columns CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and the sample name. annotation contains four columns corresponding to a four-column BED format: chromosome, gene start (0-based), gene end, and gene name.  The extent of phased variants within phased_vcf_data can restrict the extend of phased blocks returned, such that more than one phased block may be returned per annotation. Phased blocks that do not meet the thresholds of number of phased SNPs, coverage or heterzygosity as determined from input arguments are filtered out.
#' @export
#' @examples
#' data(example_phased_vcf)
#' example_phased_vcf
#' data(example_annotation)
#' example_annotation
#' createPhasedBlockList(example_phased_vcf, example_annotation)

createPhasedBlockList <- function(phased_vcf_data, annotation, min_phased_sites=2, min_coverage=10, heter=2){
  unfiltered_phased_block_list = createUnfilteredPhasedBlockList(phased_vcf_data, annotation)
  return(filterPhasedBlockList(unfiltered_phased_block_list, min_phased_sites, min_coverage, heter, phased_vcf_data))
}


#' Estimating the overdispersion parameter using a binomial normal distribution
#' @param phased_block_list A list of data.frames, each containing information for a phase block, typically a gene, as described in createPhasedBlockList
#' @return An estimation of overdispersion of the data from an idealized binomial distribution, estimated using a binomial normal distribution
#' @details The overdispersion parameter is estimated across all phased blocks using a binomial normal distribution simulated with 32-points gaussian quadrature. A value of 0 indicates the data follow the standard binomial distribution, while positive values indicate overdispersion.
#' @export
#' @examples
#' data(example_phased_block_list)
#' example_phased_block_list
#' phi.estimate(example_phased_block_list)
#' # Note: this produces an unreliable overestimate due to the small size of example_block_list

phi.estimate <- function(phased_block_list)
{
  if(all(sapply(phased_block_list,is.vector)==T)) {  # collapse allele counts within each phased block
    a=t(sapply(phased_block_list,function(q)q[5:6]))
  }else{
    a=t(sapply(phased_block_list,function(q)colSums(q[,c(5,6)])))
  }
  y=a[,1]  # number of successes
  n=rowSums(a)
  design=matrix(rep(1,nrow(a),ncol=1))

  ELIKE1<- function(param)
  {   # Exact minus log-likelihood
    beta<- 0
    sig<- param
    elike<- 0
    eta <- design%*% beta
    p<- exp(eta)/ (1+ exp(eta))
    for (ii in 1:nrow(design))
    {
      elike <- elike - log(dbinorm32(y[ii],n[ii],p[ii],sig))
    }
    return(elike)
  }

  phi=stats::nlm(ELIKE1,p=0.1,hessian=T)$estimate
  return(phi)
}

#' Calculate an exact RPASE p-value for one phased block
#' @param phased_block A single data.frame describing a phased block, typically a gene, as described in createPhasedBlockList
#' @param phi The overdispersion estimation
#' @param quiet A logical, controls if one should print out the results while calculating
#' @return An exact p-value for one phased block, typically a gene
#' @details This function implements the RPASE method for Read-backed Phasing detection of Allele-Specific Expression.  For further methodological details, see the corresponding paper (Wang et al. 2017).  Rarely, this method can require considerable time to complete when coverage is deep and variance in coverage between sites within a block (gene) is high.  In our experience this occurs in less than 0.05% of cases, and occurs most often with deeper coverage and higher SNP densities, both of which are more likely in repetitive regions.
#' @export
#' @examples
#' data(example_phased_block_list)
#' example_phased_block_list
#' RPASE(example_phased_block_list[[1]], phi=0.3)

RPASE <- function(phased_block, phi, quiet=FALSE) {
  perm = function(n, k=2) {
    grid = matrix(rep(0:n, k), n + 1, k)
    all = expand.grid(data.frame(grid))
    sums = rowSums(all)
    all[sums == n,]
  }# permutation function

  x <- as.matrix(phased_block[, c(5, 6)]) # force the structure to be matrix
  if (all(x[, 1] == x[, 2]))
    return(0.5)

  pp <- list()
  n <- nrow(x)

  if(n==1)
  {
    pp[[1]]=x
    rr <- lapply(pp, function(x) perm(sum(x)))
    if (x[1]/sum(x) == 0.5) {
      prob <- 0.5
    } else {
      if (x[1]/sum(x) > 0.5) {
        jj <- rr[[1]][(rr[[1]][, 1]/rowSums(rr[[1]])) >=x[1]/sum(x), ]
        prob=sum(mapply(dbinorm32, x = jj[, 1], n = rowSums(jj), p = 0.5,  sig = phi))
      }
      if (x[1]/sum(x) < 0.5) {
        jj <- rr[[1]][(rr[[1]][, 1]/rowSums(rr[[1]]))<= x[1]/sum(x), ]
      }
      prob <-sum(mapply(dbinorm32, x = jj[, 1], n = rowSums(jj), p = 0.5, sig = phi)) #sum of all more extremer
    }


  }else{

    # calculate the combinations that are equally or more extremer.
    yy <- x
    for (i in 1:(n - 1)) {
      pp[[i]] <- c(min(yy[, 1]), min(yy[, 2]))
      uu <- matrix(rep(pp[[i]], c(nrow(yy), nrow(yy))), ncol = 2)
      yy <- yy - uu
      if (nrow(yy) != 2) {
        tt <- matrix(c(yy[, 1][-which(yy[, 1] == 0)[1]],
                       yy[, 2][-which(yy[, 2] == 0)[1]]),
                     ncol = 2)
        yy <- tt
      }
    }#calculate the layers
    pp[[n]] <- yy[1, ]# because the last 2 levels are special, so it need to be calculated separately
    pp[[n + 1]] <- yy[2, ]
    # now pp contains all combinations that are needed to calculte probabilities

    rr <- lapply(pp, function(.x) perm(sum(.x)))# all expressions could every layer has
    qwe <- as.matrix(do.call("rbind", rr))
    nr <- unlist(lapply(rr, nrow))
    nr.length <- length(nr)
    nr.1 <- cumsum(c(0, nr[-nr.length]))
    ss <- c(seq(nr.length - 1, 1), 1)

    ee <- expand.grid(lapply(rr, function(.x) 1:nrow(.x)))
    # all possible combinations of each level, a matrix of row number of rr generated above

    kk <- NULL

    choose.wm <- function(.y) {
      # select extreme values
      nr.2 <- matrix(rep(nr.1, rep(nrow(.y), nr.length)), ncol = nr.length)
      panel <- as.matrix(.y) + nr.2
      grand.ss <- rep(ss, nrow(.y))
      grand.panel <- as.vector(t(panel))
      grand.qwe <- qwe[grand.panel, ]
      grand.mult <- grand.qwe * grand.ss
      kk.n <- matrix(0, nrow(.y), 2)
      stride <- length(ss)
      slice <- 1:stride
      for (g in 1:nrow(.y)) {
        kk.n[g,1] <- sum(grand.mult[slice, 1])
        kk.n[g,2] <- sum(grand.mult[slice, 2])
        slice <- slice + stride
      }
      return(kk.n)
    }

    kk <- choose.wm(ee)# sum of counts for 2 alleles accoring to every possible combinations

    bio.prob <- function(.y) {
      # calculate binomial normal probabilities and slice them up across the rows
      nr.2 <- matrix(rep(nr.1, rep(nrow(.y), nr.length)), ncol = nr.length)
      panel <- as.matrix(.y) + nr.2
      grand.panel <- as.vector(t(panel))
      grand.counts <- qwe[grand.panel, ]
      grand.counts[,2] <- rowSums(grand.counts)
      grand.dbinorm <- mapply(dbinorm32, x = grand.counts[, 1], n = grand.counts[, 2], p = 0.5, sig = phi)
      # calculate products of each piece
      stride <- ncol(panel)
      slice <- 1:stride
      grand.prod <- numeric(nrow(panel))
      for (g in 1:nrow(panel)) {
        grand.prod[g] <- prod(grand.dbinorm[slice])
        slice <- slice + stride
      }
      prob.res <- sum(grand.prod)
      return(prob.res)
    }# the probability (bionomial-normal) of one particular phased block under randomization

    a.2 <- sum(x[, 1])/sum(x)
    if (a.2 == 0.5) {
      prob <- 0.5
    } else {
      if (a.2 > 0.5) {
        jj <- ee[(kk[, 1]/rowSums(kk)) >= a.2, ]
      }
      if (a.2 < 0.5) {
        jj <- ee[(kk[, 1]/rowSums(kk)) <= a.2, ]
      }
      prob <- bio.prob(jj)# sum of all equal and more extremer probability from given data
    }
  }
  if (! quiet) {
      print(c(as.character(phased_block[1,4]),prob))
      utils::flush.console()
  }
  return(prob)
}


#' Calculate exact RPASE p-values for each element of a list of phased blocks
#' @param phased_block_list A list of data.frames each describing a phased block, typically a gene. This data structure is created by createPhasedBlockList
#' @param phi The overdispersion estimate; by default, it is estimated from the phased_block_list
#' @param quiet A logical, controls if one should print out the results while calculating
#' @return A vector of exact p-values
#' @details This function implements the RPASE method for detecting allele-specific expression using read-backed phasing by evaluating each phased block in phased_block_list in sequence.  Alternatively, this procedure could be parallelized by using the RPASE function itself together with functions from the parallel package or another such package.
#' @export
#' @examples
#' data(example_phased_block_list)
#' runRPASE(example_phased_block_list, phi=0.3)

runRPASE <- function(phased_block_list, phi=phi.estimate(phased_block_list), quiet=FALSE) {
  sapply(phased_block_list, RPASE, phi, quiet)
}



#!/usr/bin/Rscript --vanilla

# samples Anc/Der alleles from a Beta-Binomial model
# m = number of sites
# n = number of copies
# a,b = parameters for prior
pop_sim = function(m,n,a,b) {
    # sample from a beta-binomial
    p = rbeta(m,a,b)
    x = rbinom(m,n,p)
    y = n-x
    # construct vectors of x 1s and y 0s
    xx = lapply(x,rep,x=1)
    yy = lapply(y,rep,x=0)
    v = mapply(c,xx,yy)
    # randomize allele order per site
    v = apply(v,2,sample)
}

# header for VCF format
header =
'##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1,length=%d>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s'

vcf_str = function(data) {
    # sample ref and alt labels
    nsites = ncol(data)
    ref = replicate(nsites,sample(c("A","C","G","T"),2))
    # count AN and AC statistics
    AN = nrow(data)
    AC = apply(data==1,2,sum)

    # finalize header
    out = sprintf(header,nsites, paste0(sprintf("Indiv%05d",seq_len(nrow(data)/2)),collapse="\t"))
    
    # Generate fixed information and INFO column
    fixed = paste("1",seq_len(nsites),".",ref[1,],ref[2,],90,"PASS",sep="\t")
    info = paste0("AN=",AN,";AC=",AC,";AA=",ref[1,])

    # Generate format columns
    format = apply(data, 2, function(v) {
        m = matrix(v,ncol=2)
        m = apply(m,1,sort)
        paste0(m[1,],"/",m[2,],collapse="\t")
    })

    # Put it all together and return string
    vcf = paste(fixed,info,"GT",format,sep="\t")
    paste0(c(out,vcf),collapse="\n")
}

if(!interactive()) {
    # Parse command line
    library('getopt')
    spec = matrix(c(
        'num_sites','m',1,'integer',
        'pop_size','n',1,'integer',
        'param_a','a',1,'double',
        'param_b','b',1,'double'
    ),byrow=TRUE,ncol=4)
    opt = getopt(spec)

    # Run simulation and print output
    data = pop_sim(opt$num_sites, 2*opt$pop_size, opt$param_a, opt$param_b)
    vcf = vcf_str(data)
    cat(vcf)

    # Exit nicely
    q(status=0)
}

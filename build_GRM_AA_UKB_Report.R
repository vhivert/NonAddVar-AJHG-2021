##########################################################################
#
# Built Additive-by-Additive GRM from Additive GRM
#
# Valentin Hivert
# 04/11/2019
# Note : Partly inspired from a script of Florian Rohart
#########################################################################

##Arguments
#
# 1st arguments : number of chunks to divide the input GRM and compute the GRM AA from each part.
# 2nd argument : prefix of the input (example : toto.grm, without the .bin)
# 3rd argument : prefix of the output


temp = commandArgs(TRUE) # PART
PART = as.numeric(temp[1])
prefix = temp[2]
prefix.save = temp[3]

options(scipen=999) # to get no scientific notation (e.g. no 2e+5), which screw up the filename


Built_GRMAA = function(prefix="data_test",prefix.save="data_test_AA",size=4, power=1, PART=10){
    ##
    ## Built GRM Additive-by-Additive, following Vitezica & al. 2017 (Genetics)
    ## The GRM AA is defined as the Hadamard product of the GRMA divide by the mean diagonal of the resulting matrix
    ## Hence, we first extract the diagonal (square of the GRM A trace) and then produce the Hadamard product of the GRMA, and then divide it by the mean trace
    ## and write the final result.
    ##
    if(PART>999) stop("PART must be lower than 1000")
    sum_i <- function(i){
        return(sum(1:i))
    }
    
    ## open file connections and read in data
    BinFileName <- paste(prefix,".grm.bin",sep="") # (it is a binary file which contains the lower triangle elements of the GRM)
    IDFileName <- paste(prefix,".grm.id",sep="") #(no header line; columns are family ID and individual ID, see above).
    id <- read.table(IDFileName)
    n <- dim(id)[1] #Sample size
    
    DiagId=rep(0,n) #Vector with the position of the trace values in the GRM bin file
    for(i in 1:n){DiagId[i]=sum_i(i)}
    
    BinFile <- file(BinFileName, "rb")
    nn=n*(n+1)/2#length(grm2)
    grm_per_part = ceiling(nn/PART)
    
    
    print("extracting trace of GRM A")
    grm_current = 1
    ntrace=sum_trace2=0
    for(part in 1:PART){
        
        print(part)
        ## open file connections and read in data
        number_grm = min(grm_per_part, nn- grm_per_part*(part-1)) 
        
        cat("reading grm\n")
        grm_part <- readBin(BinFile, n=number_grm, what=numeric(0), size=size) # read only number_grm from the file
        dum=seq(grm_current,grm_current+number_grm-1)%in%DiagId
        trace=grm_part[dum==T]
        ntrace=ntrace+length(trace)
        sum_trace2=sum_trace2+sum(trace^2)
        grm_current = grm_current + number_grm
    }
    close(BinFile)
    if(ntrace!=n){ stop("The size of the trace must be equal to the number of individuals") }
    mean_trace=sum_trace2/ntrace
    rm(trace)
    rm(ntrace)
    rm(sum_trace2)
    
    BinFile <- file(BinFileName, "rb")
    print("Computing Additive-by-Additive GRM")
    for(part in 1:PART){
        
        print(part)
        ## open file connections and read in data
        if(part<10){
            bin.file.name.save <- paste(prefix.save,".part_",PART,"_00",part,".grm.bin",sep="")
        }else if (part<100){
            bin.file.name.save <- paste(prefix.save,".part_",PART,"_0",part,".grm.bin",sep="")
        }else if (part<1000){
            bin.file.name.save <- paste(prefix.save,".part_",PART,"_",part,".grm.bin",sep="")
        }
        
        number_grm = min(grm_per_part, nn- grm_per_part*(part-1)) #grm we take in that part: the max except last part where we may have to take less

        cat("reading grm\n")
        grm_part <- readBin(BinFile, n=number_grm, what=numeric(0), size=size) # read only number_grm from the file
        grm_part_power = (grm_part^power)/mean_trace
        BinFile.save <- file(bin.file.name.save, "wb")
        writeBin(con = BinFile.save, grm_part_power, size = size)
        close(BinFile.save)
        grm_current = grm_current + number_grm
    }
    close(BinFile)
    cat("done reading grm\n")

    
}

#Build GRMAA from GRMA
Built_GRMAA(prefix = prefix, prefix.save = prefix.save, power=2, PART=PART)


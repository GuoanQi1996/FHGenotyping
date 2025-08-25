# packages requirement
options(warn = -1)
necessary=c('stats','data.table')
installed=necessary %in% installed.packages()[, 'Package']
if (length(necessary[!installed]) >=1)
  install.packages(necessary[!installed],repos='http://cran.us.r-project.org')
info=lapply(necessary, function(x){library(x,character.only = T,warn.conflicts = F,quietly = T)})
options(warn = 0)

# define functions
getopt = function (spec=NULL,opt=NULL,command=get_Rscript_filename(),usage=FALSE,debug=FALSE) {
  
  # littler compatibility - map argv vector to opt
  if (is.null(opt)) {
    if (exists("argv", where = .GlobalEnv, inherits = FALSE)) {
      opt = get("argv", envir = .GlobalEnv) #nocov
    } else {
      opt = commandArgs(TRUE)
    }
  }
  
  ncol=4
  maxcol=6
  col.long.name    = 1
  col.short.name   = 2
  col.has.argument = 3
  col.mode         = 4
  col.description  = 5
  
  flag.no.argument = 0
  flag.required.argument = 1
  flag.optional.argument = 2
  
  result = list()
  result$ARGS = vector(mode="character")
  
  #no spec.  fail.
  if ( is.null(spec) ) {
    stop('argument "spec" must be non-null.')
    
    #spec is not a matrix.  attempt to coerce, if possible.  issue a warning.
  } else if ( !is.matrix(spec) ) {
    if ( length(spec)/4 == as.integer(length(spec)/4) ) {
      warning('argument "spec" was coerced to a 4-column (row-major) matrix.  use a matrix to prevent the coercion')
      spec = matrix( spec, ncol=ncol, byrow=TRUE )
    } else {
      stop('argument "spec" must be a matrix, or a character vector with length divisible by 4, rtfm.')
    }
    
    #spec is a matrix, but it has too few columns.
  } else if ( dim(spec)[2] < ncol ) {
    stop(paste('"spec" should have at least ', ncol, ' columns.',sep=''))
    
    #spec is a matrix, but it has too many columns.
  } else if ( dim(spec)[2] > maxcol ) {
    stop(paste('"spec" should have no more than ', maxcol, ' columns.',sep=''))
    
    #spec is a matrix, and it has some optional columns.
  } else if ( dim(spec)[2] != ncol ) {
    ncol = dim(spec)[2]
  }
  
  #sanity check.  make sure long names are unique, and short names are unique.
  if ( length(unique(spec[,col.long.name])) != length(spec[,col.long.name]) ) {
    stop(paste('redundant long names for flags (column ',col.long.name,' of spec matrix).',sep=''))
  }
  if ( length(stats::na.omit(unique(spec[,col.short.name]))) != length(stats::na.omit(spec[,col.short.name])) ) {
    stop(paste('redundant short names for flags (column ',col.short.name,' of spec matrix).',sep=''))
  }
  # convert numeric type to double type
  spec[,4] <- gsub("numeric", "double", spec[,4])
  
  # if usage=TRUE, don't process opt, but generate a usage string from the data in spec
  if ( usage ) {
    ret = ''
    ret = paste(ret,"Usage: ",command,sep='')
    for ( j in 1:(dim(spec))[1] ) {
      ret = paste(ret,' [-[-',spec[j,col.long.name],'|',spec[j,col.short.name],']',sep='')
      if (spec[j,col.has.argument] == flag.no.argument) {
        ret = paste(ret,']',sep='')
      } else if (spec[j,col.has.argument] == flag.required.argument) {
        ret = paste(ret,' <',spec[j,col.mode],'>]',sep='')
      } else if (spec[j,col.has.argument] == flag.optional.argument) {
        ret = paste(ret,' [<',spec[j,col.mode],'>]]',sep='')
      }
    }
    # include usage strings
    if ( ncol >= 5 ) {
      max.long = max(apply(cbind(spec[,col.long.name]),1,function(x)length(strsplit(x,'')[[1]])))
      ret = paste(ret,"\n",sep='')
      for (j in 1:(dim(spec))[1] ) {
        ret = paste(ret,sprintf(paste("    -%s|--%-",max.long,"s    %s\n",sep=''),
                                spec[j,col.short.name],spec[j,col.long.name],spec[j,col.description]
        ),sep='')
      }
    }
    else {
      ret = paste(ret,"\n",sep='')
    }
    return(ret)
  }
  
  #XXX check spec validity here.  e.g. column three should be convertible to integer
  
  i = 1
  
  while ( i <= length(opt) ) {
    if ( debug ) print(paste("processing",opt[i]))
    
    current.flag = 0 #XXX use NA
    optstring = opt[i]
    
    
    #long flag
    if ( substr(optstring, 1, 2) == '--' ) {
      if ( debug ) print(paste("  long option:",opt[i]))
      
      optstring = substring(optstring,3)
      
      this.flag = NA
      this.argument = NA
      kv = strsplit(optstring, '=')[[1]]
      # if ( !is.na(kv[2]) ) {
      if ( grepl('=', optstring) ) {
        this.flag = kv[1]
        this.argument = paste(kv[-1], collapse="=")
      } else {
        this.flag = optstring
      }
      
      rowmatch = grep( this.flag, spec[,col.long.name],fixed=TRUE )
      
      #long flag is invalid, matches no options
      if ( length(rowmatch) == 0 ) {
        stop(paste('long flag "', this.flag, '" is invalid', sep=''))
        
        #long flag is ambiguous, matches too many options
      } else if ( length(rowmatch) > 1 ) {
        # check if there is an exact match and use that
        rowmatch = which(this.flag == spec[,col.long.name])
        if(length(rowmatch) == 0) {
          stop(paste('long flag "', this.flag, '" is ambiguous', sep=''))
        }
      }
      
      #if we have an argument
      if ( !is.na(this.argument) ) {
        #if we can't accept the argument, bail out
        if ( spec[rowmatch, col.has.argument] == flag.no.argument ) {
          stop(paste('long flag "', this.flag, '" accepts no arguments', sep=''))
          
          #otherwise assign the argument to the flag
        } else {
          mode = spec[rowmatch, col.mode]
          warning_msg <- tryCatch(storage.mode(this.argument) <- mode,
                                  warning = function(w) {warning(paste(mode, "expected, got", dQuote(this.argument)))})
          if( is.na(this.argument) && !grepl("expected, got",  warning_msg) ) {
            warning(paste('long flag', this.flag, 'given a bad argument'))
          }
          result[spec[rowmatch, col.long.name]] = this.argument
          i = i + 1
          next
        }
        
        #otherwise, we don't have an argument
      } else {
        #if we require an argument, bail out
        ###if ( spec[rowmatch, col.has.argument] == flag.required.argument ) {
        ###  stop(paste('long flag "', this.flag, '" requires an argument', sep=''))
        
        #long flag has no attached argument. set flag as present.  set current.flag so we can peek ahead later and consume the argument if it's there
        ###} else {
        result[spec[rowmatch, col.long.name]] = TRUE
        current.flag = rowmatch
        ###}
      }
      
      #short flag(s)
    } else if ( substr(optstring, 1, 1) == '-' ) {
      if ( debug ) print(paste("  short option:",opt[i]))
      
      these.flags = strsplit(optstring,'')[[1]]
      
      done = FALSE
      for ( j in 2:length(these.flags) ) {
        this.flag = these.flags[j]
        rowmatch = grep( this.flag, spec[,col.short.name],fixed=TRUE )
        
        #short flag is invalid, matches no options
        if ( length(rowmatch) == 0 ) {
          stop(paste('short flag "', this.flag, '" is invalid', sep=''))
          
          #short flag has an argument, but is not the last in a compound flag string
        } else if ( j < length(these.flags) & spec[rowmatch,col.has.argument] == flag.required.argument ) {
          stop(paste('short flag "', this.flag, '" requires an argument, but has none', sep=''))
          
          #short flag has no argument, flag it as present
        } else if ( spec[rowmatch,col.has.argument] == flag.no.argument ) {
          result[spec[rowmatch, col.long.name]] = TRUE
          done = TRUE
          
          #can't definitively process this flag yet, need to see if next option is an argument or not
        } else {
          result[spec[rowmatch, col.long.name]] = TRUE
          current.flag = rowmatch
          done = FALSE
        }
      }
      if ( done ) {
        i = i + 1
        next
      }
    }
    
    #invalid opt
    if ( current.flag == 0 ) {
      stop(paste('"', optstring, '" is not a valid option, or does not support an argument', sep=''))
      #TBD support for positional args
      #if ( debug ) print(paste('"', optstring, '" not a valid option.  It is appended to getopt(...)$ARGS', sep=''))
      #result$ARGS = append(result$ARGS, optstring)
      
      # some dangling flag, handle it
    } else if ( current.flag > 0 ) {
      if ( debug ) print('    dangling flag')
      if ( length(opt) > i ) {
        peek.optstring = opt[i + 1]
        if ( debug ) print(paste('      peeking ahead at: "',peek.optstring,'"',sep=''))
        
        #got an argument.  attach it, increment the index, and move on to the next option.  we don't allow arguments beginning with '-' UNLESS
        #specfile indicates the value is an "integer" or "double", in which case we allow a leading dash (and verify trailing digits/decimals).
        if ( substr(peek.optstring, 1, 1) != '-' |
             #match negative double
             ( substr(peek.optstring, 1, 1) == '-'
               & regexpr("^-[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", peek.optstring) > 0
               & spec[current.flag, col.mode]== 'double'
             ) |
             #match negative integer
             ( substr(peek.optstring, 1, 1) == '-'
               & regexpr("^-[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", peek.optstring) > 0
               & spec[current.flag, col.mode]== 'integer'
             )
        ) {
          if ( debug ) print(paste('        consuming argument *',peek.optstring,'*',sep=''))
          
          # storage.mode(peek.optstring) = spec[current.flag, col.mode]
          mode = spec[current.flag, col.mode]
          tryCatch(storage.mode(peek.optstring) <- mode,
                   warning = function(w) {warning(paste(mode, "expected, got", dQuote(peek.optstring)))})
          result[spec[current.flag, col.long.name]] = peek.optstring
          i = i + 1
          
          #a lone dash
        } else if ( substr(peek.optstring, 1, 1) == '-' & length(strsplit(peek.optstring,'')[[1]]) == 1 ) {
          if ( debug ) print('        consuming "lone dash" argument')
          # storage.mode(peek.optstring) = spec[current.flag, col.mode]
          mode = spec[current.flag, col.mode]
          tryCatch(storage.mode(peek.optstring) <- mode,
                   warning = function(w) {warning(paste(mode, "expected, got", dQuote(peek.optstring)))}) #nocov 
          result[spec[current.flag, col.long.name]] = peek.optstring
          i = i + 1
          
          #no argument
        } else {
          if ( debug ) print('        no argument!')
          
          #if we require an argument, bail out
          if ( spec[current.flag, col.has.argument] == flag.required.argument ) {
            stop(paste('flag "', this.flag, '" requires an argument', sep=''))
            
            #otherwise set flag as present.
          } else if (
            spec[current.flag, col.has.argument] == flag.optional.argument |
            spec[current.flag, col.has.argument] == flag.no.argument 
          ) {
            x = TRUE
            storage.mode(x) = spec[current.flag, col.mode]
            result[spec[current.flag, col.long.name]] = x
          } else {
            stop(paste("This should never happen.", #nocov
                       "Is your spec argument correct?  Maybe you forgot to set", #nocov
                       "ncol=4, byrow=TRUE in your matrix call?")) #nocov
          }
        }
        #trailing flag without required argument
      } else if ( spec[current.flag, col.has.argument] == flag.required.argument ) {
        stop(paste('flag "', this.flag, '" requires an argument', sep=''))
        
        #trailing flag without optional argument
      } else if ( spec[current.flag, col.has.argument] == flag.optional.argument ) {
        x = TRUE
        storage.mode(x) = spec[current.flag, col.mode]
        result[spec[current.flag, col.long.name]] = x
        
        #trailing flag without argument
      } else if ( spec[current.flag, col.has.argument] == flag.no.argument ) {
        x = TRUE
        storage.mode(x) = spec[current.flag, col.mode]
        result[spec[current.flag, col.long.name]] = x
      } else {
        stop("this should never happen (2).  please inform the author.") #nocov
      }
    } #no dangling flag, nothing to do.
    
    i = i+1
  }
  return(result)
}
get_Rscript_filename = function() {
  prog <- sub("--file=", "", grep("--file=", commandArgs(), value=TRUE)[1])
  if( .Platform$OS.type == "windows") { 
    prog <- gsub("\\\\", "\\\\\\\\", prog)
  }
  prog
}
print_header = function(){
  version="V1.0"
  cat("###################################################\n")
  cat("## FHGenotyping (FHG)\n")
  cat(paste0("## Version: ", version,"\n"))
  cat("## Written by: Guo-An Qi, Zhejiang University\n")
  cat("## Bug report: guoan.qi@foxmail.com\n")
  cat("###################################################\n\n")
}

command = matrix(c("tfile","t",1,"character", "plink variant-major additive component genotype file",
                   "ann","a",1,"character", "Annotation for the provided variants",
                   "mis","m",1,"character","A text file specifies the type of non-synonymous variants",
                   "syn","s",1,"character","A text file specifies the type of synonymous variants",
                   "out","o",2,"character","specify prefix of output results, optional, defaulted FH",
                   "help","h",0,"logical", "parameters input instruction"),
                 byrow=T,ncol=5)
args = getopt(spec = command)

if (!is.null(args$help) || is.null(args$tfile) || is.null(args$mis) || is.null(args$syn) || is.null(args$ann)) {
  print_header()
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

# grab arguements
tfile = args$tfile
ann = args$ann
mis = args$mis
syn = args$syn

out = "FH"
if (!is.null(args$out)) {
  out = args$out
}

print_header()
cat(paste0("Options: \n  --tfile ",tfile," \n  --ann ",ann," \n  --mis ",mis," \n  --syn ",syn," \n  --out ",out,"\n\n"))


library(data.table)

GT = fread(tfile,header=T,data.table=F,nrows=20)
classes=sapply(GT, class)
classes[c(3,5,6)] = "NULL"
GT = fread(tfile,header=T,data.table=F,colClasses = classes)

GT_Annotation = fread(ann,header=T)
misTy = fread(mis,header=F)
synTy = fread(syn,header=F)

GT = merge(GT,GT_Annotation,"SNP")
colIndice = c(2,1,3,ncol(GT)-1,ncol(GT),4:(ncol(GT)-2))
GT = GT[,colIndice]
GT = GT[order(GT$Gene),]

GT.var = apply(GT[,-c(1:5)],1,var,na.rm=T)
rmIndex = which(is.na(GT.var) | GT.var==0)
if(length(rmIndex)>0){
  GT = GT[-rmIndex,]
  cat(paste0("\nWarning! \n",length(rmIndex)," variants with no polymorphisms were removed, ", nrow(GT), " variants retained.\n"))
  
}

cat(paste0(nrow(GT[which(!is.na(match(GT$Type,misTy$V1))),])," missense variants in ", length(unique(GT$Gene[which(!is.na(match(GT$Type,misTy$V1)))]))," genes remain for FH genotyping.\n"))
geneList = unique(GT$Gene)
N = ncol(GT)-5
hapNumber = c()
EH = c()
KA = c()
KS = c()

summaryDT1 = data.frame()
summaryDT2 = data.frame()
for(i in 1:length(geneList)){
  gene = geneList[i]
  dt = GT[which(GT$Gene==gene),]
  dt[is.na(dt)]=0
  missen_dt = dt[which(!is.na(match(dt$Type,misTy$V1))),]
  rows = nrow(missen_dt)
  if(rows != 0){
    if(rows > 1){
      hapmap_samplewise = apply(missen_dt[,-c(1:5)],2,paste0,collapse = "")
      hapmap = table(hapmap_samplewise)
    } else {
      hapmap_samplewise = as.character(missen_dt[,-c(1:5)])
      hapmap = table(hapmap_samplewise)
    }
    # number of haplotypes
    hpn = length(hapmap)
    hapNumber = c(hapNumber, hpn)
    # haplotypes
    haplotype_name = paste0(gene,"_",sprintf("HAP%03d", 1:hpn))
    haplotype = names(hapmap)
    haplotype_count = as.numeric(hapmap)
    
    haplotype_matrix = matrix(0,nrow = length(haplotype_name), ncol = ncol(dt)-5)
    indexRow = match(hapmap_samplewise,haplotype)
    indexCol = 1:(ncol(missen_dt)-5)
    haplotype_matrix[cbind(indexRow,indexCol)] = 1
    rownames(haplotype_matrix) = haplotype_name
    colnames(haplotype_matrix) = colnames(missen_dt)[-c(1:5)]
    hapGT = apply(haplotype_matrix,2,function(x){
      which(x==1)
    })
    summaryDT1 = rbind(summaryDT1,as.numeric(hapGT))
    
    hapSum = data.frame(Name = haplotype_name,Haplotype = haplotype, Number = haplotype_count)
    summaryDT2 = rbind(summaryDT2,hapSum)
    
    # EH
    p = hapmap/N
    eh = sum(p*log(p))*(-1)/log(N)
    EH = c(EH, eh)
  } else {
    summaryDT1 = rbind(summaryDT1,rep(1,N))
    hapNumber = c(hapNumber, 0)
    EH = c(EH, 0)
  }
  
  # Ka Ks
  missen_index = which(!is.na(match(dt$Type,misTy$V1)))
  synon_index = which(!is.na(match(dt$Type,synTy$V1)))
  if(length(missen_index)==0){
    Ka = 0
  } else {
    missense_dt = dt[missen_index,]
    hapmap_samplewise = apply(missense_dt[,-c(1:5)],2,paste0,collapse = "")
    hapmap = table(hapmap_samplewise)
    MsplitHAPs = t(do.call("cbind",strsplit(names(hapmap),split = "")))
    MsplitHAPs = matrix(apply(MsplitHAPs,2,as.integer),nrow = length(hapmap),byrow = F)
    MNi = as.numeric(hapmap)
    if(nrow(MsplitHAPs)<2){
      MHAP_SNPs = sum(MsplitHAPs!=0)
      MHAP_Sites = sum(MsplitHAPs * MNi)
    } else {
      MHAP_SNPs = rowSums(MsplitHAPs!=0)
      MHAP_Sites = rowSums(MsplitHAPs * MNi)
    }
    Ka = sum(MNi*MHAP_SNPs/MHAP_Sites,na.rm = T)/N
  }
  KA = c(KA,Ka)
  
  if(length(synon_index)==0){
    Ks = 0
  } else {
    synon_dt = dt[synon_index,]
    hapmap_samplewise = apply(synon_dt[,-c(1:5)],2,paste0,collapse = "")
    hapmap = table(hapmap_samplewise)
    SsplitHAPs = t(do.call("cbind",strsplit(names(hapmap),split = "")))
    SsplitHAPs = matrix(apply(SsplitHAPs,2,as.integer),nrow = length(hapmap),byrow = F)
    SNi = as.numeric(hapmap)
    if(nrow(SsplitHAPs)<2){
      SHAP_SNPs = sum(SsplitHAPs!=0)
      SHAP_Sites = sum(SsplitHAPs * SNi)
    } else {
      SHAP_SNPs = rowSums(SsplitHAPs!=0)
      SHAP_Sites = rowSums(SsplitHAPs * SNi)
    }
    Ks = sum(SNi*SHAP_SNPs/SHAP_Sites,na.rm = T)/N
  }
  KS = c(KS,Ks)
}
colnames(summaryDT1) = colnames(GT)[-c(1:5)]
summaryDT1 = cbind(Gene=geneList,summaryDT1)

write.table(summaryDT2,paste0(out, ".HapSTAT.txt"),sep = '\t',quote=F,col.names = T,row.names = F)
write.table(summaryDT1,paste0(out, ".GT.txt"),sep = '\t',quote=F,col.names = T,row.names = F)

geneLevelSummary = data.frame(Gene=geneList,HapN=hapNumber,EH=EH,Ka=KA,Ks=KS,`Ka/Ks`=KA/KS)
write.table(geneLevelSummary,paste0(out, ".GeneSTAT.txt"),sep = '\t',quote=F,col.names = T,row.names = F)

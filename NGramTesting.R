# NGramming Processing Script (c) by Teague Henry
#   
# NGramming Processing Script is licensed under a
# Creative Commons Attribution 3.0 Unported License.
# 
# You should have received a copy of the license along with this
# work.  If not, see <http://creativecommons.org/licenses/by/3.0/>.

#Above is a real quick license to slap on the work, use if you want, just attribute these functions
#(just the code in the document) to me

install.packages("quanteda")
install.packages("data.table")
library(quanteda)
library(data.table)

#Fair warning, this function will take some time on large corpuses. If you have a cluster I suggest you
#upload this function set to run. Cutoff processing will be much faster than p-value-testing



#TextList is a character vector of documents, each element being a document with tokens separated
#by a single white space
#Cutoff is for the threshold bigramming, only bigrams with total counts above the threshold are maintained
#P.value.test activates the p-value test, which assess p(word2|word1) NULL being random distribution of words.
#P.value.cut is the target alpha
#BF.corr activates Bonferroni correction
#dfm.return activates the return of a document frequency matrix, if False, returns processed character vector
#Return.word.list returns only the list of bigrams that are over the threshold, and does not
#                 process the data, use this to evaluate the performance of the algorithm
#Concat determines what character is used to combine bigrams
textToBigrams <- function(textList, cutoff = 100, p.value.test = F, p.value.cut = .05, 
                          dfm.return = F, return.word.list = F, concat = "_", BF.corr = F){
  corp <- corpus(textList)
  
  if(!p.value.test){
    corpDTMbigram <- dfm(corp, ngram =2, concatenator = concat)
    count <- colSums(corpDTMbigram)
    bigramsCount <-  sort(count,decreasing = T)
    if(return.word.list){
      return(bigramsCount[which(bigramsCount > cutoff)])
    }
    bigrams <- names(bigramsCount)[which(bigramsCount > cutoff)]
    if(length(bigrams) > 0){
    bigramsSpaced <- gsub(pattern = concat, replacement = " ",x = bigrams, fixed = T)

    for(i in 1:length(bigrams)){
      cat("ngram ", i, "\n")
      textList <- gsub(x = textList, pattern = paste(" ",bigramsSpaced[i], " ",sep = ""), replacement = paste(" ",bigrams[i], " ",sep = ""))
    }
    }else{
      cat("No ngrams found at this cutoff.")
    }
    
    if(dfm.return){
      return(list(dfm(corpus(textList)),bigramsCount))
    }else{
      return(list(textList, bigramsCount))
    }
  }else{
    require(data.table)
    corpDTMUnigram <- dfm(corp, ngram = 1)
    corpDTMbigram <- dfm(corp, ngram = 2, concatenator = concat)
    
    unigramCounts <-  colSums(corpDTMUnigram)
    bigramCounts <- colSums(corpDTMbigram)
    totalBigramCounts <- sum(bigramCounts)
    bigramCounts <- bigramCounts[which(bigramCounts > 1)]
    bigram.names <- do.call("rbind",strsplit(names(bigramCounts), concat, T))
    unigram.names <- names(unigramCounts)
    unigramTable <- data.table(unigram.names, unigramCounts)
    setkey(unigramTable, unigram.names)
    
    bigramTable <- data.table(bigram.names[,1], bigram.names[,2], bigramCounts, pvalue = 0)
    
    totalUnigramCount <- sum(unigramCounts)

    
    if(dim(bigramTable)[[1]] > 0){
      print(dim(bigramTable))
      for(i in 1:dim(bigramTable)[[1]]){
        prop1 <- unigramTable[.(bigramTable[i,V1]), unigramCounts]
        prop2 <- unigramTable[.(bigramTable[i,V2]), unigramCounts]/(totalUnigramCount-prop1)

        
        pval <- pbinom(bigramTable[i,bigramCounts]-1, size = prop1, prob = prop2, lower.tail = F)
        bigramTable[i, pvalue := pval]
      }
      if(return.word.list){return(bigramTable)}
      if(BF.corr){
      bigrams <- names(bigramCounts[which(bigramTable[, pvalue] < p.value.cut/dim(bigramTable)[[1]])])
      bigrams <- bigrams[order(bigramTable[which(bigramTable[, pvalue] < p.value.cut/dim(bigramTable)[[1]]),pvalue])]
      }else{  
      bigrams <- names(bigramCounts[which(bigramTable[, pvalue] < p.value.cut)])
      bigrams <- bigrams[order(bigramTable[which(bigramTable[, pvalue] < p.value.cut),pvalue])]
      }
      if(length(bigrams) > 0){
        bigramsSpaced <- gsub(pattern = concat, replacement = " ",x = bigrams, fixed = T)
        
        for(i in 1:length(bigrams)){
          cat("ngram ", i, "\n")
          textList <- gsub(x = textList, pattern = paste(" ",bigramsSpaced[i], " ",sep = ""), replacement = paste(" ",bigrams[i], " ",sep = ""))
        }
      }else{
        cat("No ngrams found at this p-value.")
      } 
      if(dfm.return){
        return(list(dfm(corpus(textList)),bigramTable))
      }else{
        return(list(textList, bigramTable))
      }
    }else{
      print("No ngram counts above 1")
    if(dfm.return){
      return(list(dfm(corpus(textList)),bigramTable))
    }else{
      return(list(textList, bigramTable))
    }}
    
  }

}


#This applies the bigrams function twice, which allows for trigrams and quadgrams (2 bigrams stuck together)
#Parameters are self explanatory

textToBigramsTrigrams <- function(textList, bigramCutoff, trigramCutoff, 
                                  p.value.test = F, p.value.cut = .05, dfm.return = F, 
                                  return.word.list  = F, BF.corr = F){
  textList <- textToBigrams(textList, bigramCutoff, p.value.test, p.value.cut, F, F, concat = "_", BF.corr = BF.corr)[[1]]
  textList <- textToBigrams(textList, trigramCutoff, p.value.test, p.value.cut, dfm.return, return.word.list, concat= ".", BF.corr = BF.corr)
  return(textList)
}

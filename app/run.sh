#!/bin/bash

#Creating filtered lists with alpha/beta chains in convenient format
#./process.sh

#Pushing information from lists above to reference classes (AlphaChain/BetaChain) removes those which dont have V(D)J regions and CDR3 and then prints to alpha_test.txt and beta_test.txt in "generated" folder
#java -classpath ./target/classes com.antigenomics.vdjstruct.Main
Rscript src/main/rscripts/cluster.R

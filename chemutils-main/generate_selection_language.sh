# Script for regenerating the python code associated with the selection language grammar

cd chemutils/molecule/selection
curl -O https://www.antlr.org/download/antlr-4.13.2-complete.jar
java -jar antlr-4.13.2-complete.jar -Dlanguage=Python3 -visitor SelectionLanguage.g4 -o generated
rm antlr-4.13.2-complete.jar
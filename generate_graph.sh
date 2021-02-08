# Script to generate a call graph for Mikado using pyan3

pip install pyan3

cd Mikado/
cat <(echo __init__.py) <(find _transcripts configuration daijin daijin loci parsers picking preparation scales serializers subprograms transcripts utilities  -iname "*py") | xargs pyan3 --dot --no-defines -G -a --colored -e -g > ../Mikado.dot
cd ../
cat Mikado.dot | dot -Tsvg -Granksep=1.5 > Mikado.svg

#!/bin/bash

APCOCOA_PATH="/home/garak/Dokumente/Uni/Promotion/ApCoCoA/apcocoa.org/ApCoCoALib/trunk"

mkdir tmp-build;
cd tmp-build;
echo -n > include
echo -n > namespace
echo -n > c-file
echo -n > h-file
FILES="CoeffField TMonoid TOrdering Term Polynomial F4Algorithm F4Utils";
for j in ${FILES}; do
	i="../../src/${j}.C";
	grep "^[[:space:]]*#include[[:space:]]*<" $i >> include;
	grep "namespace" $i >> namespace;
	grep -v "#include" $i | grep -v "using namespace" >> c-file;
done;

for j in ${FILES}; do
	i="../../include/${j}.H";
	grep "^[[:space:]]*#include[[:space:]]*<" $i >> include;
	grep "namespace" $i >> namespace;
	grep -v "#include" $i | grep -v "using namespace" >> h-file;
done;

sort include | uniq >> header.H
sort namespace | uniq >> header.H
#cat h-file >> ParallelF4.H
#cat ../footer.H >> ParallelF4.H

sed -e "/<<<INCLUDE>>>/r header.H" -e "/<<<INCLUDE>>>/d" ../target.H > tmp.H
sed -e "/<<<HEADER>>>/r h-file" -e "/<<<HEADER>>>/d" tmp.H > ParallelF4.H
sed -e "/<<<SOURCE>>>/r c-file" -e "/<<<SOURCE>>>/d" ../target.C > ParallelF4.C

rm include namespace c-file h-file header.H tmp.H
#mv ParallelF4.C ParallelF4.H ..
mv ParallelF4.C "${APCOCOA_PATH}/src/Algebraic/"
mv ParallelF4.H "${APCOCOA_PATH}/include/ApCoCoA/Algebraic/"
cd ..
rmdir tmp-build

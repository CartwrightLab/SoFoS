#!/bin/bash
#
# library(rmutil)
# f = function(n,a,b) dbetabinom(0:n,n,a/(a+b),a+b)

echo "Running Integration Test: unfolded+tagged ..."

SOFOS="${1:-../sofos}"
VCF=test-default.vcf

output=$(${SOFOS} -q -a 0.1 -b 1 -n 10 -u -t ${VCF})

hnum=$(echo "${output}" | grep '^#' | wc -l)
if [[ 10 -ne "${hnum}" ]]; then
	echo "  ERROR: Heading has ${hnum} lines instead of 10."
	exit 1
fi
header=$(echo "${output}" | awk 'NR==11')
if [[ $header != "Number,Prior,Observed,Posterior" ]]; then
	echo "  ERROR: Header has wrong column names."
	exit 1
fi

# Check posterior
expected="0.649602923876743
0.476864186047936
0.323085024970344
0.207787973205851
0.128175439263581
0.0772737404468384
0.0475690799391035
0.0320738639692642
0.0246226484295746
0.0198606972724853
0.0130844225782785"

echo "${output}" | grep -v '^#' | awk -F, 'NR>1 {print $4}' \
	| paste -d, - <(echo "${expected}") \
	| awk -F, -v m=1e-7 'function abs(v) {return v < 0 ? -v : v} {a=abs(($1-$2))/$2; if(a > m){m = a} } END{exit (m > 1e-7) }' 
if [[ $? -ne 0 ]]; then
	echo "  ERROR: Maximum relative error in posterior exceeds 1e-7."
	exit 1
fi

# Check prior
expected="1.50324356453683
0.150324356453683
0.0826783960495252
0.0578748772346676
0.0448530298568674
0.0367794844826313
0.0312625618102366
0.0272430895774919
0.0241782420000241
0.0217604178000216
0.0198019801980198"

echo "${output}" | grep -v '^#' | awk -F, 'NR>1 {print $2}' \
	| paste -d, - <(echo "${expected}") \
	| awk -F, -v m=1e-7 'function abs(v) {return v < 0 ? -v : v} {a=abs(($1-$2))/$2; if(a > m){m = a} } END{exit (m > 1e-7) }' 
if [[ $? -ne 0 ]]; then
	echo "  ERROR: Maximum relative error in prior exceeds 1e-7."
	exit 1
fi

# Check obs
expected="0
0
1.9
0
0
0
0
0
0.1
0
0"

echo "${output}" | grep -v '^#' | awk -F, 'NR>1 {print $3}' \
	| paste -d, - <(echo "${expected}") \
	| awk -F, -v m=1e-7 'function abs(v) {return v < 0 ? -v : v} {a=abs(($1-$2)); if(a > m){m = a} } END{exit (m > 1e-7) }' 
if [[ $? -ne 0 ]]; then
	echo "  ERROR: Maximum error in obs exceeds 1e-7."
	exit 1
fi

# Check indexes
expected="0
1
2
3
4
5
6
7
8
9
10"

echo "${output}" | grep -v '^#' | awk -F, 'NR>1 {print $1}' \
	| paste -d, - <(echo "${expected}") \
	| awk -F, -v m=1e-7 'function abs(v) {return v < 0 ? -v : v} {a=abs(($1-$2)); if(a > m){m = a} } END{exit (m > 1e-7) }' 
if [[ $? -ne 0 ]]; then
	echo "  ERROR: Maximum error in obs exceeds 1e-7."
	exit 1
fi

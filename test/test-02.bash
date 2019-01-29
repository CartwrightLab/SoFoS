#!/bin/bash
#
# library(rmutil)
# f = function(n,a,b) dbetabinom(0:n,n,a/(a+b),a+b)

echo "Running Integration Test: folded+refalt ..."

SOFOS="${1:-../sofos}"
VCF=test-default.vcf

output=$(${SOFOS} -q -a 0.1 -b 1 -n 10 -f -r ${VCF})

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

nrow=$(echo "${output}" | grep -v '^#' | awk -F, 'END {print NR}')
if [[ 7 -ne "${nrow}" ]]; then
	echo "  ERROR: Body has wrong number of rows."
	exit 1
fi

# Check posterior
expected="0.613042662838015
0.481731622579739
0.354864195501706
0.259140144535515
0.200677862510758
0.0905435120342663"

echo "${output}" | grep -v '^#' | awk -F, 'NR>1 {print $4}' \
	| paste -d, - <(echo "${expected}") \
	| awk -F, -v m=1e-7 'function abs(v) {return v < 0 ? -v : v} {a=abs(($1-$2))/$2; if(a > m){m = a} } END{exit (m > 1e-7) }' 
if [[ $? -ne 0 ]]; then
	echo "  ERROR: Maximum relative error in posterior exceeds 1e-7."
	exit 1
fi

# Check prior
expected="1.52304554473485
0.172084774253705
0.106856638049549
0.0851179668121595
0.076115591667104
0.0367794844826313"

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
2
0
0
0
"

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
5"

echo "${output}" | grep -v '^#' | awk -F, 'NR>1 {print $1}' \
	| paste -d, - <(echo "${expected}") \
	| awk -F, -v m=1e-7 'function abs(v) {return v < 0 ? -v : v} {a=abs(($1-$2)); if(a > m){m = a} } END{exit (m > 1e-7) }' 
if [[ $? -ne 0 ]]; then
	echo "  ERROR: Maximum error in obs exceeds 1e-7."
	exit 1
fi

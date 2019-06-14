echo 'First 10 components of theta_star'
grep star $1 | cut -d' ' -f4| cut -d',' -f1-10

echo ''

echo 'Top 10 presentations'

grep game $1 | cut -d' ' -f3 \
    | awk '{
    split($0, arr, ",");
    asort(arr);
    for (i=1; i<=length(arr); i++) {
        printf "%s,", arr[i] };
        printf RS
    }' | sort -n | uniq -c | sort -rn | head -n 10


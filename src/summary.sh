grep largest $1 | cut -d' ' -f5 \
    | awk '{
    split($0, arr, ",");
    asort(arr);
    for (i=1; i<=length(arr); i++) {
        printf "%s,", arr[i] };
        printf RS
    }' | sort -n | uniq -c | sort -rn | head -n $2


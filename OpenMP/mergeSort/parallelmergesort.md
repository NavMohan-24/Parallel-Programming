## Psuedocode for Parallel Merge Sort

```
/**
 * A: Input array
 * B: Output array
 * lo: lower bound
 * hi: upper bound
 * off: offset
 **/


algorithm parallelMergesort(A, lo, hi, B, off) is
    len := hi - lo + 1
    if len == 1 then
        B[off] := A[lo]
    else let T[1..len] be a new array
        mid := ⌊(lo + hi) / 2⌋ 
        mid' := mid - lo + 1
        fork parallelMergesort(A, lo, mid, T, 1)
        parallelMergesort(A, mid + 1, hi, T, mid' + 1) 
        join 
        parallelMerge(T, 1, mid', mid' + 1, len, B, off)
```
#### Psuedocode for Parallel Merge

```
algorithm merge(A[i...j], B[k...ℓ], C[p...q]) is
    inputs A, B, C : array
           i, j, k, ℓ, p, q : indices

    let m = j - i,
        n = ℓ - k

    if m < n then
        swap A and B  // ensure that A is the larger array: i, j still belong to A; k, ℓ to B
        swap m and n

    if m ≤ 0 then
        return  // base case, nothing to merge

    let r = ⌊(i + j)/2⌋
    let s = binary-search(A[r], B[k...ℓ])
    let t = p + (r - i) + (s - k)
    C[t] = A[r]

    in parallel do
        merge(A[i...r], B[k...s], C[p...t])
        merge(A[r+1...j], B[s...ℓ], C[t+1...q])
```


**References**

1. [Parallel Merge Sort](https://en.wikipedia.org/wiki/Merge_sort#Parallel_merge_sort)
2. [Parallel Merge](https://en.wikipedia.org/wiki/Merge_algorithm#Parallel_merge)
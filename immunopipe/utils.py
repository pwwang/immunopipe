def chunk_list(array, lens):
    start = 0
    out = []
    for length in lens:
        out.append(array[start: (start + length)])
        start += length
    return out

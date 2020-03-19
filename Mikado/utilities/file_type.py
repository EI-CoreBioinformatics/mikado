import sys

magic_dict = {
    b"\x1f\x8b\x08": b"application/gzip",
    b"\x42\x5a\x68": b"application/x-bzip2",
    b"\x50\x4b\x03\x04": b"application/zip"
    }


max_len = max(len(x) for x in magic_dict)


def filetype(filename):
    with open(filename, "rb") as f:
        file_start = f.read(max_len)
    for magic, ftype in magic_dict.items():
        if file_start.startswith(magic):
            return ftype
    return b"application/txt"

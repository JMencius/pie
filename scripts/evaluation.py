import ctypes
import os
from scripts.myblock import myblock

##N50, NG50, NG90, Switch Error, Generalize hamming distance, phasing percentage



# load C so

script_dir = os.path.dirname(os.path.realpath(__file__))
lib = ctypes.CDLL(f"{script_dir}/hamming.so")

# hamming distance C function
lib.hamming_distance.restype = ctypes.c_int
lib.hamming_distance.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
def c_hamming_distance(str1, str2) -> int:
    if len(str1) != len(str2):
        raise ValueError("Length of the two input strings should be equal.")
    return lib.hamming_distance(str1.encode(), str2.encode())



# hamming comparision C function
lib.hamming_comparison.restype = ctypes.POINTER(ctypes.c_int)
lib.hamming_comparison.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
# free result C function
lib.free_result.argtypes = [ctypes.POINTER(ctypes.c_int)]
lib.free_result.restype = None
def hamming_comparison(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Length of the two input strings should be equal.")
    
    length = len(str1)
    
    result_ptr = lib.hamming_comparison(str1.encode(), str2.encode(), length)
    result = [result_ptr[i] for i in range(length)]
    lib.free_result(result_ptr)
    
    return result



def blockwise_evaluate(in_block) -> tuple:
    query: str = ''.join(in_block.left)
    truth: str = ''.join(in_block.truthleft)
    print(hamming_comparison(query, truth))
    print(c_hamming_distance(query, truth))








"""
if __name__ == "__main__":
    str1 = "1@001"
    str2 = "10011"
    distance = c_hamming_distance(str1, str2)
    print(f"Hamming distance: {distance}")
"""



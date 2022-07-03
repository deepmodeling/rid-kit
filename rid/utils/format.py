from typing import List

def list_to_string(
        input_list: List, 
        split_sign: str
    ) -> str:
    input_list = [str(element) for element in input_list]
    return split_sign.join(input_list)

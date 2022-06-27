from typing import List

def list_to_string(
        input_list: List, 
        split_sign: str
    ) -> str:
    return split_sign.join(input_list)

from typing import Optional, Tuple

dict1 = {"a": 1, "b": 2}
genome_list = [genome for genome in dict1]
print(genome_list)
print(next(iter(dict1))[0])
for item in dict1.items():
    print(item)

def func(**kwargs):
    param1 = kwargs.get("param1")
    param2 = kwargs.get("param2")
    param3 = kwargs.get("param3")
    print(param1, param2, param3)

func(param1=1, param2=2)

second: Optional[Tuple[str, int]] = None
print(second[1])

# def validate_type(param: Any, param_type: Type, message: str = None) -> None:
#     """
#     This function checks if the given parameter is of the given type
#     :param param: parameter to be checked
#     :param param_type: type of the parameter
#     :param message: message to be printed if the parameter is not of the given type
#     :return: None if the parameter is of the given type, Raises TypeError otherwise
#     """
#     if not isinstance(param, param_type):
#         message = message if message is not None else "Parameter must be of type {}".format(
#             param_type)
#         raise TypeError(message)


# with open(
#         "C:\\Users\\omers\\Documents\\CS\Year1\Semester1\Intro2CS\Intro2cs_Psuedo_Aligner\json\\ref_mid_clean_k31_dump.json",
#         "r") as eg_file:
#     example_json = json.load(eg_file)
#     reference_json = json.loads(reference.to_json())
#     print(example_json == reference_json)
#
#
# def write_to_kdb_file(kdb_file_path: str,
#                       kmer_reference: KmerReference) -> bool:
#     """
#     This function handles the writing of the kmer_reference object to our own kdb file
#     The kmer_reference will be pickled and its content stored in a gzip file
#     :param kdb_file_path: the path to the kdb file
#     :param kmer_reference: the kmer reference object to write
#     :return: True if the kdb file was written, False otherwise
#     """
#     try:
#         validate_file_type(kdb_file_path, program_constants.KDB_FILE_TYPES)
#     except TypeError as e:
#         print(e)
#         return False
#     pickled_object = pickle.dumps(kmer_reference)
#     try:
#         with open(kdb_file_path, 'wb') as kdb_file:
#             kdb_file.write(gzip.compress(pickled_object))
#     except (FileNotFoundError, PermissionError, IOError) as e:
#         print("An error occurred while writing to the given reference file.")
#         return False
#     except Exception as e:
#         print(
#             "An unexpected error occurred while writing to the given reference file.")
#         return False
#     return True
#
#
# def decompress_kdb_file(kdb_file_path: str) -> Optional[KmerReference]:
#     """
#     This function receives a path to a kdb file, checks its integrity and decompresses the kdb file
#     as well as unpickles it
#     :param kdb_file_path: the path to the kdb file
#     :return: the KmerReference within the kdb_file
#     """
#     try:
#         validate_file_type(kdb_file_path, program_constants.KDB_FILE_TYPES)
#     except TypeError as e:
#         print(e)
#         return None
#     try:
#         with open(kdb_file_path, 'rb') as kdb_file:
#             pickled_object = gzip.decompress(kdb_file.read())
#     except (FileNotFoundError, PermissionError, IOError) as e:
#         print("An error occurred while reading the given kdb file.")
#         return None
#     except (ValueError, OSError) as e:
#         print("There was an error while decompressing the kdb file")
#         return None
#     except Exception as e:
#         print(
#             "An unexpected error occurred while reading the given kdb file.")
#         return None
#
#     try:
#         kmer_reference = pickle.loads(pickled_object)
#     except (UnpicklingError, TypeError) as e:
#         print(
#             "There was a problem unpickling the kdb file decompressed content")
#         print(e)
#         return None
#     except Exception as e:
#         print(
#             "There was an unexpected problem unpickling the kdb file decompressed content")
#         return None
#     return kmer_reference
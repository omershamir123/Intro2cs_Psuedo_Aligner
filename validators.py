# FILE : validators.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: a validator module - containing various validation methods
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from typing import Any, Type, Iterable, List


def validate_not_empty(param: Any) -> bool:
    """
    This function checks whether a parameter equals to None
    :param param: the parameter to be checked
    :return: None if not None, otherwise raises ValueError
    """
    return param is not None

def validate_file_type(file_path: str, file_types: List[str], message: str = None) -> None:
    """
    This function receives a file path and makes sure it's of the same type
    :param file_path: string containing the file path
    :param file_types: a list containing the file types that are valid for the file
    :param message: message to be printed if the file is not of the given file_type
    :return: None if the file is in the correct type, otherwise raises a TypeError
    """
    if not validate_not_empty(file_path):
        raise TypeError("The parameter file_path is empty")
    if not any(file_path.endswith(file_type) for file_type in file_types):
        message = message if message is not None else "File {} must be of types {}".format(
            file_path, file_types)
        raise TypeError(message)


def validate_type(param: Any, param_type: Type, message: str = None) -> None:
    """
    This function checks if the given parameter is of the given type
    :param param: parameter to be checked
    :param param_type: type of the parameter
    :param message: message to be printed if the parameter is not of the given type
    :return: None if the parameter is of the given type, Raises TypeError otherwise
    """
    if not isinstance(param, param_type):
        message = message if message is not None else "Parameter must be of type {}".format(
            param_type)
        raise TypeError(message)


def validate_above_value(param: Any, threshold: Any, allow_equality: bool,
                         message: str = None) -> None:
    """
    This function checks if the given parameter is above the given threshold
    :param param: parameter to be checked
    :param threshold: threshold to be checked
    :param allow_equality: can the parameter be equal to or not to the threshold
    :param message: message to be printed if the parameter is not of the given type
    :return: None if the parameter is of the given type, Raises ValueError otherwise
    """
    if not validate_not_empty(param):
        raise TypeError("The parameter is empty to be checked that it's above {} is None".format(threshold))
    if not allow_equality:
        if param <= threshold:
            message = message if message is not None else "Parameter must be greater than {}".format(
                threshold)
            raise ValueError(message)
    else:
        if param < threshold:
            message = message if message is not None else "Parameter must be greater than or equal to {}".format(
                threshold)
            raise ValueError(message)


def validate_values_in_given_list(iter: Iterable, allowed_values: List[Any], message: str = None) -> None:
    """
    This function checks if the given items in an iterable are all part of the allowed_values
    :param iter: the iterable to go over
    :param allowed_values: the allowed values
    :param message: the message to be printed in case of an error
    :return: None if the items in an iterable are all part of the allowed_values, Raises ValueError otherwise
    """
    valid = all(item in allowed_values for item in iter)
    if not valid:
        message = message if message is not None else "All items must be in {}".format(allowed_values)
        raise ValueError(message)


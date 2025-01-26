# FILE : validators.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: a validator module - containing various validation methods
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from typing import Any, Type, Iterable, List


def general_file_validation(file_path: str, actions_needed: str) -> None:
    """
    This function receives a file path and makes sure it's a valid file path
    that can be opened and read
    :param file_path: string containing the file path
    :param actions_needed: a string containing "artwb" being the possible actions when opening a file
    :return: None if the file is readable, Raises ValueError otherwise
    """
    pass


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

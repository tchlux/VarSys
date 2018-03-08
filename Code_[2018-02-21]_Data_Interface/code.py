# Function takes a python list "data" as input and uses advanced
# statistical techniques to determine if more experiments need to be
# run in order for a sample to accurately capture the behavior of the
# system at the given configuration.
def should_stop(current_data):
    if len(current_data) >= 50:
        return True
    else:
        return False
    # ^^ We're working on something more advanced, but this will be
    # the framework of what we give you.
    

# ============================
#      Checking CSV Files     
# ============================

CONFIG_SEPARATOR = "--" # Separates the name of the parameter and its values
CONFIG_DELIMITER = "," # Separates the individual values for a parameter
DEFAULT_REQUIRE_ALL_COLUMNS = True # True if all columns must be provided
# Custom exceptions for specific errors
class BadConfigFile(Exception):       pass
class BadColumnInData(Exception):     pass
class MissingColumnInData(Exception): pass
class BadRowInData(Exception):        pass
class BadValueInData(Exception):      pass


# Given a path to a CSV file and a path to an allowed configurations
# file, check to make sure the contents of a provided CSV are
# valid. If the contents are invalid, raise a descriptive error.
def check_csv(csv_path, config_path, 
                 require_all_columns=DEFAULT_REQUIRE_ALL_COLUMNS):
    # Open the configuration file and read the allowed configurations
    with open(config_path) as conf_file:
        allowed_columns = [line.strip() for line in conf_file.readlines()]
        column_names_and_values = {}
        for column_specification in allowed_columns:
            # Check to make sure the separator exists in the the row
            if CONFIG_SEPARATOR not in column_specification:
                raise(BadConfigFile(("\n\n'%s' is missing the required "+
                                     "separator string '%s'.")%(
                                         column_specification, CONFIG_SEPARATOR)))
            seperator_index = column_specification.index(CONFIG_SEPARATOR)
            # Get the name of the column (before the config separator)
            column_name = column_specification[:seperator_index].strip()
            # Get the values of the column (comma separated)
            column_values = column_specification[
                seperator_index+len(CONFIG_SEPARATOR):]
            column_values = [val.strip() for val in column_values.split(CONFIG_DELIMITER)]
            # Store this in a dictionary for lookup on validation
            column_names_and_values[column_name] = column_values

    # Open the data file and validate the data values
    with open(csv_path) as data_file:
        all_lines = data_file.readlines()
        # Retrieve the header to know the order of columns
        header = [col_name.strip() for col_name in all_lines.pop(0).split(CONFIG_DELIMITER)]
        # Verify that each header is a valid header.
        for col_name in header:
            if col_name not in column_names_and_values:
                raise(BadColumnInData((
                    "\n\nThe column name '%s' is not specified as an allowed"+
                    " configuration parameter by '%s'.")%(col_name, config_path)))
        # Check that all columns were provided, if that's expected
        if require_all_columns:
            for col_name in column_names_and_values:
                if col_name not in header:
                    raise(MissingColumnInData((
                        "\n\nThe column name '%s' is specified by '%s' but "+
                        "is not provided in '%s'.")%(
                            col_name, config_path, csv_path)))
        # Process each row of data
        for i, line in enumerate(all_lines):
            line = [col_val.strip() for col_val in line.strip().split(CONFIG_DELIMITER)]
            # Check that the line has the correct number of rows.
            if len(line) != len(header):
                raise(BadRowInData(("\n\nWrong number of columns on row %i in '%s'."+
                                    " Expected %i, found %i.")%(
                                       i+2, csv_path, len(header), len(line))))
            # Check that all values are acceptable
            for (col_val, col_name) in zip(line, header):
                if col_val not in column_names_and_values[col_name]:
                    raise(BadValueInData(("\n\nRow %i in '%s' has value '%s' for column"+
                                          " '%s', which is not a member of %s.")%(
                                              i+2, csv_path, col_val, col_name, 
                                              column_names_and_values[col_name])))
        # The CSV must be a valid CSV.
        return True

if __name__ == "__main__":
    # Test the function on our sample CSV file
    check_csv("input-sample.csv", "sample-types.config")
    check_csv("bad-sample.csv", "sample-types.config")

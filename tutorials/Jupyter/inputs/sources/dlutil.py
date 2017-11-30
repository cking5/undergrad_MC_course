"""Utility methods related to DL-MONTE"""

import os

def remove_comments(line):

    """Remove comments (#) from FIELD, CONFIG, CONTROL lines

    Also compress any mulitple white spaces.
    """

    line = line.split("#")[0]
    line = " ".join(line.split())
    return line.strip()


def load_ascii(infile, directory=os.curdir):

    """Load ASCII input file and return lines as list of strings

    Read and remove comments and blank lines from an ASCII inout
    file (FIELD, CONFIG, CONTROL).

    Returns:
      A list of lines with actual content

    Raises:
      IOError: there's a problem with the file
    """

    filename = os.path.join(directory, infile)

    # Read and store (non-comment, non-blank) lines
    # The order of the lines is important to be able to
    # parse the input correctly

    filehandle = open(filename, "r")

    with filehandle:
        lines = filehandle.readlines()

    content = []

    for line in lines:
        line = line.strip()
        line = remove_comments(line)
        if line != "":
            content.append(line)

    return content

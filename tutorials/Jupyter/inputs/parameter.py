"""Utility to store int or float  parameters with a label

The label is added merely for information: there is no
extra functionality associated with it.

The module has a factory class Parameter(value, label)
which returns either ParameterInt ot ParameterFloat
depending on type(value).
"""

class ParameterInt(int):

    """Integer parameter with a label"""

    def __new__(cls, *args):

        """Arguements:
        args[0] is the value
        """

        return super(ParameterInt, cls).__new__(cls, args[0])


    def __init__(self, value, label):

        """Arguments
        value (integer):              the value
        label (string, Label, ...):   associated label
        """

        super(ParameterInt, self).__init__()
        self.label = label


    def __repr__(self):

        """Return string including the label"""

        parameter = "value= {!s}, label= {!r}".format(self, self.label)

        return "ParameterInt({!s})".format(parameter)


class ParameterFloat(float):

    """Real parameter with label"""

    def __new__(cls, *args):

        """Arguments:
        args[0] is the value
        """

        return super(ParameterFloat, cls).__new__(cls, args[0])


    def __init__(self, value, label):

        """Arguments
        value (float):                the value
        label (string, Label, ...):   a description
        """

        super(ParameterFloat, self).__init__()
        self.label = label


    def __repr__(self):

        """Return a string including the label"""

        parameter = "value= {!s}, label= {!r}".format(self, self.label)

        return "ParameterFloat({!s})".format(parameter)


class Parameter(object):

    """A factory to return a Parameter of the correct type"""

    def __new__(cls, value, label):

        """Arguments:
        value (int or float)
        label (string or Label)
        """

        if isinstance(value, int):
            return ParameterInt(value, label)
        elif isinstance(value, float):
            return ParameterFloat(value, label)
        else:
            raise TypeError("A Parameter is either int or float")

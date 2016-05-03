#! /usr/bin/python


DEFAULT_PARAMS = {'--version': 0
}

DEFAULT_DELIMITERS = ' \n\t'

# TODO: Allow definition of parameter types, and then parse them according to types
# TODO: Allow definition of function that would process each parameter (or all of them)
# TODO: Create custom exceptions


# This utility class represents program parameters
# A user can define which parameter is followed by how many arguments, and then parse input
# either from a string or a list of comand line arguments.
# dictionary params contains all allowed parameters annd a number of arguments for each
class Parser:
    def __init__(self, params = DEFAULT_PARAMS):
        self.params = params


    def parseString(self, string, delimiters = DEFAULT_DELIMITERS):
        args = string.split(deimiters)
        return self.parseCmdArgs(args)


    def parseCmdArgs(self, args):
        paramsdict = {}
        i = 0
        while i < len(args):
            param = args[i]
            if param not in self.params:
                raise Exception('Invalid parameter %s' % param)

            paramlist = []
            for k in xrange(self.params[param]):
                i += 1
                if i >= len(args):
                    raise Exception('Parameter %s does\'t have enough arguments (%d/%d)' % (aparamrg, k, self.params[param]))
                paramlist.append(args[i])
            paramsdict[param] = paramlist
            i += 1

        return paramsdict

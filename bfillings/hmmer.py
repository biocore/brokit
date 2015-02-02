# -----------------------------------------------------------------------------
# Copyright (c) 2015--, biocore development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

""" Application controller for HMMER v. 3.1b1 """

from os.path import abspath, join, dirname

from burrito.parameters import ValuedParameter, FlagParameter
from burrito.util import (CommandLineApplication, ResultPath,
                          ApplicationError)


class hmmer(CommandLineApplication):
    """ HMMER generic application controller.

    	Instead of instantiating this class, instantiate a
    	subclass for each binary. Available subclasses are:

    	hmmer_nhmmer: DNA homology search with profile HMMs
    """

    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameters(self, data):
        """ Set the input path (a fasta filepath)
        """
        # The list of values which can be passed on a per-run basis
        allowed_values = ['-o', '--noali', '--tblout',
                          '--acc', '-E', '--cpu']

        unsupported_parameters = set(data.keys()) - set(allowed_values)
        if unsupported_parameters:
            raise ApplicationError(
                "Unsupported parameter(s) passed when calling nhmmer: %s" %
                ' '.join(unsupported_parameters))

        for v in allowed_values:
            # turn the parameter off so subsequent runs are not
            # affected by parameter settings from previous runs
            self.Parameters[v].off()
            if v in data:
                # turn the parameter on if specified by the user
                self.Parameters[v].on(data[v])

        return ''

    def getHelp(self):
	    """ Method that points to documentation """
	    help_str = """
	    HMMER package is hosted at:
	    http://hmmer.janelia.org/software

	    The following papers should be cited if this resource is used:

	    [1] Finn R.D., Clements J. and Eddy S.R., "HMMER web server:
	        interactive sequence similarity searching", Nucl Acids Res,
	        2011, Web Server Issue 39:W29-W37.
	    """
	    return help_str


class hmmer_nhmmer(hmmer):
	""" Controls the nhmmer command of HMMER.
		nhmmer implements DNA homology search with profile HMMs.
	"""

    _command_ = 'nhmmer'
    _parameters_ = {
        # Output results to specified file
        '-o': ValuedParameter('-', Name='o', Delimiter=' ',
                              IsPath=True),
        # Don't output alignments
        '--noali': FlagParameter('--', Name='noali', Value=True),
        # Save parseable table of hits to file
        '--tblout': ValuedParameter('--', Name='tblout', Delimiter=' ',
                                    IsPath=True),
        # Prefer accessions over names in output
        '--acc': FlagParameter('--', Name='acc', Value=True),
        # Report sequences <= this E-value threshold in output
        '-E': ValuedParameter('-', Name='E', Delimiter=' ',
                              IsPath=False, Value=1e-5),
        # Number of parallel CPU working to use for multithreads
        '--cpu': ValuedParameter('--', Name='cpu', Delimiter=' ',
        		                 IsPath=False, Value=1)
    }

    def _get_result_paths(self, data):
        """ Set the result paths """

        result = {}

        # Direct output to file
        result['FastaOutput'] = ResultPath(
            Path=self.Parameters['-o'].Value,
            IsWritten=self.Parameters['-o'].isOn())

        # Parseable table of hits
        result['TableOutput'] = ResultPath(
            Path=self.Parameters['--tblout'].Value,
            IsWritten=self.Parameters['--tblout'].isOn())

    def getHelp(self):
        """ Method that points to documentation """
        help_str = """
        nhmmer (newly part of HMMER 3.1b1 package) is hosted at:
        http://hmmer.janelia.org/software

        The following papers should be cited if this resource is used:

        [1] Wheeler T.J. and Eddy S.R., "nhmmer: DNA homology search
            with profile HMMs", Bioinformatics, 2013, 29(19):2487-9.
        """
        return help_str


def run_nhmmer(fasta_filepath,
			   profile_hmm,
			   output_filepath=None,):
	""" Search a DNA model or alignment against a
	    DNA database

	    Parameters
	    ----------

	    Return
	    ------
	"""


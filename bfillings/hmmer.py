# -----------------------------------------------------------------------------
# Copyright (c) 2015--, biocore development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

""" Application controller for HMMER v. 3.1b1 """

from skbio.parse.sequences import parse_fasta

from burrito.parameters import ValuedParameter, FlagParameter
from burrito.util import (CommandLineApplication, ResultPath)


class hmmer(CommandLineApplication):
    """ HMMER generic application controller.

        Instead of instantiating this class, instantiate a
        subclass for each binary. Available subclasses are:

        hmmer_nhmmer: DNA homology search with profile HMMs
    """

    _suppress_stdout = False
    _suppress_stderr = False

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

    _command = 'nhmmer'
    _input_handler = '_input_as_string'
    _parameters = {
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
        result['AlignmentsOutput'] = ResultPath(
            Path=self.Parameters['-o'].Value,
            IsWritten=self.Parameters['-o'].isOn())

        # Parseable table of hits
        result['TableOutput'] = ResultPath(
            Path=self.Parameters['--tblout'].Value,
            IsWritten=self.Parameters['--tblout'].isOn())

        return result

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
               output_filepath,
               tblout_filepath=None,
               report_artifacts=False,
               noali=True,
               acc=True,
               evalue=None,
               threads=1,
               HALT_EXEC=False):
    """ Search a DNA model or alignment against a
        DNA database

        Parameters
        ----------
        fasta_filepath : string
            filepath to FASTA sequences to search against
            profile HMM
        profile_hmm : string
            filepath to HMM profile
        output_filepath : string
            filepath to output alignments
        tblout_filepath : string, optional
            filepath to parseable table of hits
        report_artifacts : boolean, optional
            if true, return a list of artifact labels
        noali : boolean, optional
            if false, do not output alignments
        acc : boolean, optional
            if true, prefer accessions over names in output
        evalue : float, optional
            E-value threshold for reporting matching sequences
        threads : integer, optional
            number of threads to use

        Return
        ------
        app_result : dict
            dictionary of open file objects
        artifacts : list
            list containing artifact labels (can be empty)
    """
    # Instantiate the object
    nhmmer = hmmer_nhmmer(HALT_EXEC=HALT_EXEC)

    # Set parameters
    if not noali:
        nhmmer.Parameters['--noali'].off()
    if not acc:
        nhmmer.Parameters['--acc'].off()
    if evalue:
        nhmmer.Parameters['-E'].on(evalue)
    if threads > 0:
        nhmmer.Parameters['--cpu'].on(threads)
    else:
        raise ValueError("Number of threads must be a positive integer.")
    if tblout_filepath:
        nhmmer.Parameters['--tblout'].on(tblout_filepath)

    # artifacts can be deduced using the default output
    # (set by '-o') but the file is more difficult to parse
    # than '--tblout' and contains more information than
    # we need for this step
    if (report_artifacts and tblout_filepath is None):
        raise ValueError("A table of parseable table of hits is required "
                         "as output in order to deduce the artifacts.")

    nhmmer.Parameters['-o'].on(output_filepath)

    data = "%s %s" % (profile_hmm, fasta_filepath)
    app_result = nhmmer(data)

    artifacts = []

    # Find all artifacts
    if report_artifacts:
        # Get list of all hits
        hits = set()
        for line in app_result['TableOutput']:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            hits.add(line[0])

        # Get list of all input sequences
        all_seqs = set()
        with open(fasta_filepath, 'U') as fasta_f:
            for label, seq in parse_fasta(fasta_f):
                all_seqs.add(label)

        # Get sequences (artifacts) which didn't match to the HMM
        artifacts = list(all_seqs - hits)

    return app_result, artifacts

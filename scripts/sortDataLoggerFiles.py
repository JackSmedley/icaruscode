#!/usr/bin/env python
#
# Documentation can be rendered in python interactive session with
# `import sortDataLoggerFiles; help(sortDataLoggerFiles)` and on command line
# with `sortDataLoggerFiles --help`.
# 
# It may require ROOT.
# 
# Changes:
# 20210411 (petrillo@slac.fnal.gov) [1.0]
#   first public version
#

import sys, os
import re
import logging

__doc__ = """Sorts a list of data logger output files.

File paths are read from all the specified file lists in sequence, or from
standard input if no file list is specified.

If a line is encountered that does not match the typical file name pattern,
that line is ignored and a warning is printed.

Comments and empty lines at the beginning of the first file list are printed
at the top of the output as they are. All other comments and empty lines are
printed at the end of the output.

Note that it is possible to sort "in place" by specifying the same file list as
input and output.

Duplicate files are files on the same run, data logger cycle and data logger
number. By default, only the first of the duplicate files is written into the
output list, and only the number of duplicates is printed. Options allow to
write a detailed list of duplicate files on screen and on disk, or not to check
for duplication altogether.

"""

__author__ = 'Gianluca Petrillo (petrillo@slac.stanford.edu)'
__date__ = 'February 26, 2021'
__version__ = '1.0'


class CycleCompareClass:
  """Provides less() to compare numbers with a single offset cycle.
  
  For example, with offset 3 the order of [0:20] would be, from the lowest:
  [ 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 0, 1, 2, ]
  """
  def __init__(self, first): self.first = first
  def less(self, a, b): return (a < b) == ((a < self.first) == (b < self.first))
# class CycleCompareClass

class FileInfoClass:
  """This class collects information about a input file, including a sorting
  criterium.
  """
  
  XRootDprotocolHead = 'root://fndca1.fnal.gov:1094/'
  XRootDprotocolDir = 'pnfs/fnal.gov/usr'
  POSIXprotocolHead = '/'
  POSIXprotocolDir = 'pnfs'
  
  Pattern = re.compile(r"data_dl(\d+)_run(\d+)_(\d+)_(.*)\.root")
  POSIXPattern = re.compile(
    POSIXprotocolHead.replace('.', r'\.')
    + POSIXprotocolDir.replace('.', r'\.')
    + r"/([^/]+)/(.*)")
  XRootDPattern = re.compile(
    XRootDprotocolHead.replace('.', r'\.')
    + XRootDprotocolDir.replace('.', r'\.')
    + r"/(.*)"
    )
  _DataLoggerSorter = CycleCompareClass(first=4)
  
  @staticmethod
  def getFirstDataLogger(index): return FileInfoClass._DataLoggerSorter.first
  @staticmethod
  def setFirstDataLogger(index):
    FileInfoClass._DataLoggerSorter = CycleCompareClass(first=index)
  
  def __init__(self,
               line: "input file line (should include endline)",
               source: "an arbitrary identifier to track the origin of the line" = None,
               ):
    """Constructor: use and parse the specified input file line."""
    self.line = line
    self.source = source
    self.path = line.strip()
    self.protocolAndDir, self.name = os.path.split(self.path)
    match = FileInfoClass.Pattern.match(self.name)
    self.is_file = match is not None
    if self.is_file: self.dataLogger, self.run, self.pass_ = map(int, match.group(1, 2, 3))
  # __init__()
  
  def __lt__(self, other):
    """Comparison: run, then pass, then (offset cycled) data logger number."""
    if not self.is_file:
      raise RuntimeError \
        ("Sorting not supported for non-file objects ('%s')" % self.path)
    # if
    if self.run < other.run: return True
    if self.run > other.run: return False
    
    if self.pass_ < other.pass_: return True
    if self.pass_ > other.pass_: return False
    
    return \
      FileInfoClass._DataLoggerSorter.less(self.dataLogger, other.dataLogger)
  # __lt__()
  
  def pathToXRootD(self) -> "stored file path in XRootD format":
    if not self.is_file:
      raise RuntimeError(
       "XRootD conversion not supported for non-file objects ('%s')" % self.path
       )
    # if not file
    match = FileInfoClass.POSIXPattern.match(self.path)
    return os.path.join(
      FileInfoClass.XRootDprotocolHead, FileInfoClass.XRootDprotocolDir,
      *match.group(1, 2)
      ) if match else self.path
    
  # pathToXRootD()
  
  def pathToPOSIX(self) -> "stored file path in POSIX (PNFS local) format":
    if not self.is_file:
      raise RuntimeError(
       "XRootD conversion not supported for non-file objects ('%s')" % self.path
       )
    # if not file
    match = FileInfoClass.XRootDPattern.match(self.path)
    return os.path.join(
        FileInfoClass.POSIXprotocolHead, FileInfoClass.POSIXprotocolDir,
        match.group(1)
      ) if match else self.path
    
  # pathToXRootD()
  
# class FileInfoClass


class MinimumAccumulator:
  def add(self, data, key = None):
    if key is None: key = data
    try:
      if key >= self.minKey: return False
    except AttributeError: pass # no self.minKey yet?
    self.minKey = key
    self.minData = data
    return True
  # add()
  def min(self): return self.minData
# class MinimumAccumulator


def findFirstCycle(files):
  firstLogger = None
  firstPassFiles = []
  wrapped = False
  for info in fileInfo:
    if firstLogger == info.dataLogger: break # cycle completed
    if wrapped and info.dataLogger > firstLogger: break # cycle completed
    
    if firstLogger is None: firstLogger = info.dataLogger
    elif not wrapped and info.dataLogger < firstLogger: wrapped = True
    
    firstPassFiles.append(info)
    logging.debug("Added cycle %d logger %d to first cycle list",
      info.pass_, info.dataLogger)
  # for
  return firstPassFiles
# findFirstCycle()


def extractFirstEvent(filePath):
  try: import ROOT
  except ImportError:
    raise RuntimeError("""ROOT python module could not be loaded.
      In this condition, you'll have to skip the autodetection of the first logger
      by explicitly specifying its number as option to the script."""
      )
  # try ... except
  logging.debug("Opening '%s' for event number check...", filePath)
  srcFile = ROOT.TFile.Open(filePath, "READ")
  if not srcFile:
    raise RuntimeError \
      ("Failed to open '%s' for event number extraction." % filePath)
  #
  firstEvent = next(iter(srcFile.Events)) # go PyROOT
  firstEventNumber = firstEvent.EventAuxiliary.event() # keep going PyROOT
  
  logging.debug("First event from '%s': %d", filePath, firstEventNumber)
  return firstEventNumber
# extractFirstEvent()


def detectFirstLogger(fileInfo):
  lowestEvent = MinimumAccumulator()
  for info in fileInfo:
    firstEvent = extractFirstEvent(info.pathToXRootD())
    lowestEvent.add(info, key=firstEvent)
  firstLogger = lowestEvent.min().dataLogger
  logging.debug("Detected first logger: %d", firstLogger)
  return firstLogger
# detectFirstLogger()


def buildFileIndex(
  fileInfo: "list with information from all files",
  ) -> "a dictionary: { key -> list of files }":
  
  fileKey = lambda info: ( info.run, info.pass_, info.dataLogger )
  index = {}
  for info in fileInfo:
    index.setdefault(fileKey(info), []).append(info)
  return index
# buildFileIndex()


if __name__ == "__main__":
  
  logging.basicConfig(level=logging.INFO)
  
  import argparse
  
  parser = argparse.ArgumentParser(description=__doc__)
  parser.set_defaults(skipDuplicates=True)
  
  parser.add_argument('inputFiles', nargs="*", metavar='inputFileNames',
    help='input file lists [one from stdin by default]')
  parser.add_argument('--firstlogger', type=int,
    help='index of the first data logger in the cycle')
  parser.add_argument('--output', '-o', default=None,
    help=
     'name of the file to write the resulting list into (overwritten!) [stdout]'
    )
  parser.add_argument('--nooutput', action="store_true",
    help='do not print on screen nor write to file the files in input')
  
  duplGroup = parser.add_argument_group(title="duplicate file options")
  duplGroup.add_argument('--printduplicates', '-d', action="store_true",
    help='print duplicate files on screen')
  duplGroup.add_argument('--skipduplicates', '-S', dest='skipDuplicates',
    action="store_true",
    help='do not include duplicate files in the list (default)'
    )
  duplGroup.add_argument('--keepduplicates', '-K', dest='skipDuplicates',
    action="store_false",
    help='include also duplicate files in the list (default)'
    )
  duplGroup.add_argument('--duplicatelist', '-D', type=str, default=None,
    help='name of a file list to be created with duplicate entries')
  
  parser.add_argument('--xrootd', '--root', '-X', action="store_true",
    help='convert the paths to XRootD URL')
  parser.add_argument('--posix', '-P', action="store_true",
    help='convert the paths to local POSIX path')
  parser.add_argument('--debug', action="store_true",
    help='prints out debugging messages')
  parser.add_argument \
    ('--version', '-V', action='version', version='%(prog)s ' + __version__)

  args = parser.parse_args()
  
  if args.debug: logging.getLogger().setLevel(logging.DEBUG)
  
  if args.xrootd and args.posix:
    raise RuntimeError("XRootD and POSIX output format options are exclusive.")
  
  printDuplicates = args.printduplicates
  skipDuplicates = args.skipDuplicates
  makeDuplicateList = args.duplicatelist
  
  # "sources" are given directly as input (None = sys.stdin)
  sources = args.inputFiles if args.inputFiles else [ "<stdin>" ]
  
  # "inputFiles" are all the files found in the sources
  inputFiles = (
      [ file_ ] if file_.endswith('.root') else open(file_, 'r')
        for file_ in args.inputFiles 
    ) if args.inputFiles else [ sys.stdin, ]
  
  # example: /pnfs/icarus/persistent/users/ascarpel/trigger/4989/decoded/17247391_0/data_dl2_run4989_1_20210219T015125_20210219T200434-decode.root
  
  preComments = []
  postComments = []
  fileInfo = []
  sourceNames = []
  for iSource, file_ in enumerate(inputFiles):
    isSingleFile = isinstance(file_, list) and len(file_) <= 1
    for iLine, line in enumerate(file_):
      info = FileInfoClass(line, source=( iSource, None if isSingleFile else iLine ))
      if not info.is_file:
        if not info.path or info.path.startswith('#'):
          (postComments if fileInfo else preComments).append(info.line)
          continue
        else:
          logging.warning \
           ("Line %d ('%s') does not match file pattern." % (iLine, info.path))
          continue
      # if not file
      fileInfo.append(info)
    # for line in file
  # for input files
  
  if fileInfo and (args.firstlogger is None):
    # uses internal FileInfoClass ordering (firstLogger not set: any will do)
    fileInfo.sort()
    firstPassFiles = findFirstCycle(fileInfo)
    assert firstPassFiles
    firstLogger = detectFirstLogger(firstPassFiles)
    
  else: firstLogger = args.firstlogger if args.firstlogger is None else 4
  
  FileInfoClass.setFirstDataLogger(firstLogger)
  
  fileInfo.sort() # uses internal FileInfoClass ordering
  
  #
  # deal with duplicates
  #
  if printDuplicates or makeDuplicateList or skipDuplicates:
    nDuplicates = 0
    fileIndex = buildFileIndex(fileInfo)
    uniqueFiles = [] if skipDuplicates else None
    duplicateFiles = [] if makeDuplicateList else None
    # we rely on insertion-ordered dictionary guarantee of Python 3.7
    for fileList in fileIndex.values():
      if uniqueFiles is not None: uniqueFiles.append(fileList[0])
      if len(fileList) > 1:
        nDuplicates += len(fileList) - 1
        if duplicateFiles is not None: duplicateFiles.extend(fileList[1:])
        if printDuplicates:
          firstSource = fileList[0].source[0]
          msg = f"Run {info.run} cycle {info.pass_} data logger {info.dataLogger} with {len(fileList) - 1} duplicates of"
          
          if len(sources) > 1: msg += f" {sources[fileList[0].source[0]]}"
          if fileList[0].source[1] is not None: msg += f" line {fileList[0].source[1]}"
          msg += ":"
          for info in fileList[1:]:
            if info.source[0] != firstSource: msg += f"{sources[info.source[0]]}"
            if info.source[1] is not None: msg += f" line {info.source[1]}"
            msg += ";"
          # for
          logging.info(msg)
        # if print duplicates
      # if duplicates
    # for
    if nDuplicates: logging.info(f"Found {nDuplicates} duplicate files.")
    if duplicateFiles:
      with open(makeDuplicateList, 'w') as DuplicateListFile:
        for info in duplicateFiles: # lines still have their <CR>
          print(info.line, file=DuplicateListFile, end='')
      logging.info(f"{nDuplicates} duplicate file names written in '{makeDuplicateList}'.")
    # if we have duplicates and we write them
  # if print or store duplicates
  
  fileListContent = uniqueFiles if skipDuplicates else fileInfo
  
  
  #
  # print everything
  #
  
  # NOTE: keep this after all the input has been read,
  #       so that input files can be safely overwritten
  if not args.nooutput:
    outputFile = open(args.output, 'w') if args.output else sys.stdout
    
    # <CR> were not removed from `line`
    for line in preComments: outputFile.write(line)
    for info in fileListContent:
      if args.posix: line = info.pathToPOSIX() + '\n'
      elif args.xrootd: line = info.pathToXRootD() + '\n'
      else: line = info.line
      outputFile.write(line)
    for line in postComments: outputFile.write(line)
    
    if outputFile is not sys.stdout:
      logging.info \
        (f"{len(fileListContent)} file entries written into '{outputFile.name}'.")
      del outputFile
    # if
  else:
    logging.info(f"Found {len(fileListContent)} file entries.")
  # if ... else
  
  sys.exit(0)
  
# main

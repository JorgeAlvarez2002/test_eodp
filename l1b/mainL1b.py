
# MAIN FUNCTION TO CALL THE L1B MODULE

from l1b.src.l1b import l1b

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = r'C:\Users\jorge\PycharmProjects\EOPD\EODP_CODE\test_eodp\auxiliary'
indir = r"C:\\Users\jorge\PycharmProjects\EOPD\EODP-TS-E2E\jorgeoutputs_END2END"
outdir = r"C:\\Users\jorge\PycharmProjects\EOPD\EODP-TS-E2E\l1b_jorgeoutputs"

# Initialise the ISM
myL1b = l1b(auxdir, indir, outdir)
myL1b.processModule()

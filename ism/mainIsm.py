
# MAIN FUNCTION TO CALL THE ISM MODULE

from ism.src.ism import ism

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = r'C:\Users\jorge\PycharmProjects\EOPD\EODP_CODE\test_eodp\auxiliary'
#parte final E2E
indir = r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-E2E\sgm_out"
outdir = r"C:\\Users\jorge\PycharmProjects\EOPD\EODP-TS-E2E\jorgeoutputs_END2END"
#para apartado OPTICAL STAGE!!!!
#-indir = r"C:\Users\jorge\PycharmProjects\EOPD\EODP-TS-ISM\input\gradient_alt100_act150"
#-outdir = r"C:\\Users\jorge\PycharmProjects\EOPD\EODP-TS-ISM\jorgeoutputs"

# Initialise the ISM
myIsm = ism(auxdir, indir, outdir)
myIsm.processModule()

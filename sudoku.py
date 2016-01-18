import argparse
import csp
import time

###############################################################
# NOTE: Don't edit this file
###############################################################
def main(args):
    ###
    # args.input: this will give you the path to game.txt file
    ###

    ##########################################################
    # backtracking
    ##########################################################
	tic = time.clock()
	solution, consistencyChecks = csp.backtracking(args.input)
	toc = time.clock()
	timeItr = toc - tic

	print "Backtracking:"
	print "Execution Time: " + str(timeItr)
	print "Consistency Checks: " + str(consistencyChecks)
	print "Solution: " + str(solution) + '\n'

    ##########################################################
    # backtracking + MRV
    ##########################################################
	tic = time.clock()
	solution, consistencyChecks = csp.backtrackingMRV(args.input)
	toc = time.clock()
	timeItr = toc - tic

	print "backtrackingMRV:"
	print "Execution Time: " + str(timeItr)
	print "Consistency Checks: " + str(consistencyChecks)
	print "Solution: " + str(solution) + '\n'

    ##########################################################
    # backtracking + MRV + fwd
    ##########################################################
	tic = time.clock()
	solution, consistencyChecks = csp.backtrackingMRVfwd(args.input)
	toc = time.clock()
	timeItr = toc - tic

	print "backtrackingMRVfwd:"
	print "Execution Time: " + str(timeItr)
	print "Consistency Checks: " + str(consistencyChecks)
	print "Solution: " + str(solution) + '\n'

    ##########################################################
    # backtracking + MRV + CP
    ##########################################################
	tic = time.clock()
	solution, consistencyChecks = csp.backtrackingMRVcp(args.input)
	toc = time.clock()
	timeItr = toc - tic

	print "backtrackingMRVcp:"
	print "Execution Time: " + str(timeItr)
	print "Consistency Checks: " + str(consistencyChecks)
	print "Solution: " + str(solution) + '\n'

    ##########################################################
    # minConflict
    ##########################################################
	tic = time.clock()
	solution, consistencyChecks = csp.minConflict(args.input)
	toc = time.clock()
	timeItr = toc - tic

	print "minConflict:"
	print "Execution Time: " + str(timeItr)
	print "Consistency Checks: " + str(consistencyChecks)
	print "Solution: " + str(solution) + '\n'


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="HomeWork Three")
	parser.add_argument("--input", type=str)
	args = parser.parse_args()
	main(args)
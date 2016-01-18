import copy
from Queue import PriorityQueue
from Queue import Queue

###########################################
# you need to implement five funcitons here
###########################################

# Board configuration
N = 0
M = 0
K = 0
goal = []
queue = PriorityQueue()
unAssignedCells = []
backTrackconsistencyChecks = 0
backTrackMRVconsistencyChecks = 0
backTrackMRVFCconsistencyChecks = 0
backTrackMRVCPconsistencyChecks = 0
minConflictConsistencyChecks = 0
board = []
BackTrackingMRVfwdcount = 0
BackTrackingMRVcount = 0
minConflictCount = 0
sudokuQueue = None
consistencyChecks = 0

def backtracking(filename):
    ###
    # use backtracking to solve sudoku puzzle here,
    # return the solution in the form of list of 
    # list as describe in the PDF with # of consistency
    # checks done
    ###
    global goal, consistencyChecks, board
    consistencyChecks = 0
    readInput(filename)
    doBackTracking(copy.deepcopy(board))
    gBoard = [[(int)for j in range(N)] for i in range(N)]
    if len(goal) != 0:
        for x in range(N):
            for y in range(N):
                gBoard[x][y] = goal[x][y].value
    else:
        print "No solution exists"
    return ([gBoard], consistencyChecks)

'''
->Here we read the board and create a 2D Matrix which contains cells(in the form of objects).
->Each object has rownumber,columnnumber,value,domain,etc field associated with each of the cell).
->Then this method calls doBackTracking. We iterate through each of the cells in the matrix and then pull out the possible domain(or values that this cell can take.
Default values are 1 to N).
->For each of these values, we check if a collision happens. If there is no collision, we go ahead with the next cell recursively.
->If something doesn't workout, we start backtracking and we choose the next possible value from the cell's domain and again go ahead with the next cell recursively.
->If all the values in the domain are exhausted with return with a false which means that we couldn't find a solution.
'''
def doBackTracking(game):
    global N, goal, consistencyChecks
    isSolved = checkGoalState(game)
    if isSolved:
        return True
    cell = None
    for x in range(0, N):
        for y in range(0, N):
            if game[x][y].value == -1:
                cell = game[x][y]
                cellDomain = cell.domain
                for i in range(0, N):
                    value = cellDomain[i]
                    isValidValueExists = False
                    if not checkCollision(game, cell, x, y, value):
                        consistencyChecks += 1
                        isValidValueExists = True
                        cell.value = value
                        if doBackTracking(game):
                            return True
                        else:
                            isValidValueExists = False
                    if i == N - 1 and not isValidValueExists:
                        cell.value = -1
                        # consistencyChecks -= 1
                        return False
    return False


def checkGoalState(game):
    for x in range(0, N):
        for y in range(0, N):
            if (game[x][y].value == -1):
                return False
    global goal
    goal = copy.deepcopy(game)
    return True

# Checks if there is constraint violation
def checkCollision(game, cell, x, y, value):
    global N, M, K
    for i in range(0, N):
        if game[i][y].id != cell.id:
            if game[i][y].value == value:
                return True
    for i in range(0, N):
        if game[x][i].id != cell.id:
            if game[x][i].value == value:
                return True
    x1 = (x / M) * M
    y1 = (y / K) * K
    for i in range(x1, x1 + M):
        for j in range(y1, y1 + K):
            if (i >= 0 and i < N and j >= 0 and j < N):
                if game[i][j].id != cell.id:
                    if game[i][j].value == value:
                        return True
    return False

# prints the board
def printBoard(game):
    for x in range(0, N):
        for y in range(0, N):
            print game[x][y].value

'''
->Here we read the board and create a 2D Matrix which contains cells(in the form of objects).
->Each object has rownumber,columnnumber,value,domain,etc field associated with each of the cell).
->Then this method calls doBacktrackingMRV. We update the updateDomains and call a makeHeuristic function which is a priorityQueue implementation.
-> In the priorityQueue we are keeping the cells and their domain length values, so that we can pull out the cell which has the least length.
->We use that cell and then we iterate through each of the cells in the matrix and then pull out the possible domain(or values that this cell can take. Default values are 1 to N).
->For each of these values, we check if a collision happens. If there is no collision, we go ahead with the next cell recursively.Again domains are updated and makeHeuristic is called.
->If something doesn't workout, we start backtracking and we choose the next possible value from the cell's domain and again go ahead with the next cell recursively just like above.
->If all the values in the domain are exhausted with return with a false which means that we couldn't find a solution.
'''
def backtrackingMRV(filename):
    ###
    # use backtracking + MRV to solve sudoku puzzle here,
    # return the solution in the form of list of 
    # list as describe in the PDF with # of consistency
    # checks done
    ###
    global goal, queue, N, board, backTrackMRVconsistencyChecks
    backTrackMRVconsistencyChecks = 0
    readInput(filename)
    queue = PriorityQueue()
    doBacktrackingMRV(copy.deepcopy(board))
    gBoard = [[(int)for j in range(N)] for i in range(N)]
    if len(goal) != 0:
        for x in range(N):
            for y in range(N):
                gBoard[x][y] = goal[x][y].value
    else:
        print "No solution exists"
    return ([gBoard], backTrackMRVconsistencyChecks)

def doBacktrackingMRV(game):
    global N, M, K, queue, backTrackMRVconsistencyChecks
    if checkGoalState(game):
        return True
    updateDomains(game)
    makeHeuristic(game)
    cell = queue.get()[1]
    domain = copy.deepcopy(cell.domain)
    for i in domain:
        if not checkCollision(game, cell, cell.row, cell.column, i):
            backTrackMRVconsistencyChecks += 1
            cell.value = i
            if checkGoalState(game):
                return True
            if doBacktrackingMRV(game):
                return True
            # backTrackMRVconsistencyChecks -= 1
            cell.value = -1;
    return False

# updates the domains of variables
def updateDomains(board):
    # global game

    actualDomain = []
    for d in range(1, N + 1):
        actualDomain.append(d)

    for i in range(0, N):
        for j in range(0, N):
            if not board[i][j].readOnly:
                board[i][j].domain = copy.deepcopy(actualDomain)

    for i in range(0, N):
        for j in range(0, N):
            if board[i][j].value == -1:
                continue
            updateDomain1(board[i][j], board)


def updateDomain1(cell, game):
    global M, K, N
    x = cell.row
    y = cell.column
    for i in range(0, N):
        if game[i][y].id != cell.id:
            if not game[i][y].readOnly:
                if cell.value in game[i][y].domain:
                    game[i][y].domain.remove(cell.value)

    for i in range(0, N):
        if game[x][i].id != cell.id:
            if not game[x][i].readOnly:
                if cell.value in game[x][i].domain:
                    game[x][i].domain.remove(cell.value)

    x1 = (x / M) * M
    y1 = (y / K) * K
    for i in range(x1, x1 + M):
        for j in range(y1, y1 + K):
            if (i >= 0 and i < N and j >= 0 and j < N):
                if game[i][j].id != cell.id:
                    if not game[i][j].readOnly:
                        if cell.value in game[i][j].domain:
                            game[i][j].domain.remove(cell.value)

def updateDomain(cell, game):
    global M, K, N
    celldomain = cell.domain
    cellRow = cell.row
    cellColumn = cell.column
    if (cell.value == -1 and celldomain != None):
        # print celldomain
        for i in range(0, N):
            if (game[i][cellColumn].value != -1):
                if int(game[i][cellColumn].value) in celldomain:
                    celldomain.remove(int(game[i][cellColumn].value))
                    # print celldomain
        for j in range(0, N):
            if (game[cellRow][j].value != -1):
                if int(game[cellRow][j].value) in celldomain:
                    celldomain.remove(int(game[cellRow][j].value))

        x1 = (cellRow / M) * M
        y1 = (cellColumn / K) * K
        for i in range(x1, x1 + M):
            for j in range(y1, y1 + K):
                if (i >= 0 and i < N and j >= 0 and j < N):
                    #if game[i][j].id != cell.id:
                    if (int(game[i][j].value) in celldomain):
                        celldomain.remove(int(game[i][j].value))

# Builds the heuristic to pick a variable with minimum remaining value
def makeHeuristic(game):
    global N, queue
    queue = PriorityQueue()
    for x in range(0, N):
        for y in range(0, N):
            if game[x][y].value != -1:
                continue
            cell = game[x][y]
            queue.put((len(cell.domain), cell))

'''
->Same as above but here after we do the collision check, we also check if the insertion of a value in the cell will result in any of the neighbouring cell's domain becoming empty.[checkifAnyDomainIsEmpty method]
->In that case,we start backtracking from there itself.
'''
def backtrackingMRVfwd(filename):
    ###
    # use backtracking +MRV + forward propogation
    # to solve sudoku puzzle here,
    # return the solution in the form of list of 
    # list as describe in the PDF with # of consistency
    # checks done
    ###
    global board, goal, queue,BackTrackingMRVfwdcount
    readInput(filename)
    queue = PriorityQueue()
    #updateDomains(board)
    game = copy.deepcopy(board)
    #makeHeuristic(game)
    doBacktrackingMRVfwd(game)
    gBoard = [[(int)for j in range(N)] for i in range(N)]
    if len(goal) != 0:
        for x in range(N):
            for y in range(N):
                gBoard[x][y] = goal[x][y].value
    else:
        print "No solution exists"
    return ([gBoard], BackTrackingMRVfwdcount)

def doBacktrackingMRVfwd(game):
    global N, M, K, queue, BackTrackingMRVfwdcount
    if checkGoalState(game):
        return True
    updateDomains(game)
    makeHeuristic(game)
    cell = queue.get()[1]
    domain = copy.deepcopy(cell.domain)
    for i in domain:
        if not checkCollision(game, cell, cell.row, cell.column, i):
            cell.value = i
            updateDomains(game)
            someDomainEmpty = checkifAnyDomainIsEmpty(game)
            if(not someDomainEmpty):
                #print "here!@Varun"
                BackTrackingMRVfwdcount += 1
                if checkGoalState(game):
                    return True
                if doBacktrackingMRVfwd(game):
                    return True
            cell.value = -1;
    return False

def checkifAnyDomainIsEmpty(game):
     for i in range(0,N):
         for j in range(0,N):
             if((game[i][j].value == -1) and (len(game[i][j].domain)==0)):
                 return True
     return False

# Implements AC3 algorithm
def backtrackingMRVcp(filename):
    ###
    # use backtracking + MRV + cp to solve sudoku puzzle here,
    # return the solution in the form of list of 
    # list as describe in the PDF with # of consistency
    # checks done
    ###
    global goal, queue, N, board, backTrackMRVCPconsistencyChecks
    backTrackMRVCPconsistencyChecks = 0
    readInput(filename)
    queue = PriorityQueue()
    getUnassignedCells(board)
    doBacktrackingMRVcp(copy.deepcopy(board))
    gBoard = [[(int)for j in range(N)] for i in range(N)]
    if len(goal) != 0:
        for x in range(N):
            for y in range(N):
                gBoard[x][y] = goal[x][y].value
    else:
        print "No solution exists"
    return (gBoard, backTrackMRVCPconsistencyChecks)


def doBacktrackingMRVcp(game):
    global N, M, K, queue, unAssignedCells, backTrackMRVCPconsistencyChecks
    if checkGoalState(game):
        return True
    updateDomains(game)
    makeHeuristic(game)
    cell = queue.get()[1]
    for i in cell.domain:
        backTrackMRVCPconsistencyChecks += 1
        domain = copy.deepcopy(cell.domain)
        propagateConstraints(unAssignedCells)
        game[cell.row][cell.column].value = i
        for c in unAssignedCells:
            if cell.id == c.id:
                unAssignedCells.remove(c)

        if doBacktrackingMRVcp(game):
           return True

        # backTrackMRVCPconsistencyChecks -= 1
        game[cell.row][cell.column].value = -1
        updateDomains(game)
        getUnassignedCells(game)

def propagateConstraints(cells):
    # Queue to maintain in arcs
    arcList = []
    for cell in cells:
        affectingCells = getAffectingCells(cell)
        for aCell in affectingCells:
            arcList.append((cell, aCell))

    for x in arcList:
        arcAB = x[0]
        arcBA = x[1]
        if maintainArcConsistency(arcAB, arcBA):
            for arc in affectingCells:
                arcList.append((arc, arcAB))
                affectingCells = getAffectingCells(arcBA)
                affectingCells.remove(arcAB)

def maintainArcConsistency(self, arcAB, arcBA):
    global unAssignedCells
    aDomain = unAssignedCells[arcAB].domain
    x = aDomain[0]
    if x in aDomain:
        bDomain = unAssignedCells[arcBA].domain
        if x in bDomain:
            bDomain.remove(x)
        if len(bDomain) == 0:
            unAssignedCells[arcAB].remove(x)
            return True
    return False

def getUnassignedCells(game):
    global unAssignedCells
    unAssignedCells[:] = []
    for i in range(N):
        for j in range(N):
            if game[i][j] == -1:
                unAssignedCells.append(game[i][j])

def getAffectingCells(cell, game):
    cells = []
    global M, K, N
    x = cell.row
    y = cell.column
    for i in range(0, N):
        if game[i][y].id != cell.id:
            if not game[i][y].readOnly:
                if game[i][y].vale == -1:
                    if game[i][y] not in cells:
                        cells.append(game[i][y])

    for i in range(0, N):
        if game[x][i].id != cell.id:
            if not game[x][i].readOnly:
                if game[i][y].vale == -1:
                    if game[x][i] not in cells:
                        cells.append(game[x][i])

    x1 = (x / M) * M
    y1 = (y / K) * K
    for i in range(x1, x1 + M):
        for j in range(y1, y1 + K):
            if (i >= 0 and i < N and j >= 0 and j < N):
                if game[i][j].id != cell.id:
                    if not game[i][j].readOnly:
                        if game[i][y].vale == -1:
                            if game[i][j] not in cells:
                                cells.append(game[i][j])

'''
->The approach is similar to the 2nd Algorithm. But, the choice of the value from the domain depends on the minConflict heuristic.
->We choose the value among the cell's domain which results in least number of conflicts due to its insertion.
->Rest of the steps are similar.
'''
def minConflict(filename):
    ###
    # use minConflict to solve sudoku puzzle here,
    # return the solution in the form of list of 
    # list as describe in the PDF with # of consistency
    # checks done
    ###
    global board, goal, queue, minConflictCount
    minConflictCount = 0
    readInput(filename)
    queue = PriorityQueue()
    updateDomains(board)
    game = copy.deepcopy(board)
    prepPriorityQueue(game)
    Successvalue = dominConflict(game)
    print Successvalue
    gBoard = [[(int)for j in range(N)] for i in range(N)]
    if len(goal) != 0:
        for x in range(N):
            for y in range(N):
                gBoard[x][y] = goal[x][y].value
    else:
        print "No solution exists"
    return ([gBoard], minConflictCount)

def dominConflict(game):
    global N, M, K, queue, minConflictCount
    if checkGoalState(game):
        return True
    updateDomains(game)
    cell = None
    for i in range(0,N):
        for j in range(0,N):
            if((game[i][j].value ==-1) and (game[i][j].domain != None)):
                cell = game[i][j]
                for dom in cell.domain:
                    value = getminConflictValue(game,cell)
                    if not checkCollision(game, cell, cell.row, cell.column, value):
                        cell.value = value
                        updateDomains(game)
                        minConflictCount += 1
                        if checkGoalState(game):
                            return True
                        if dominConflict(game):
                            return True
                        cell.value = -1;
                    updateDomains(game)
    return False

def deleteDomain(cell,game,value):
    global M,K,N
    #celldomain = cell.domain
    cellRow = cell.row
    cellColumn = cell.column
    #if(cell.value==-1 and celldomain!=None):
        #print celldomain
    for i in range(0,N):
        if((game[i][cellColumn].value ==-1)  and game[i][cellColumn].domain !=None ):
            cellRowDomain = game[i][cellColumn].domain
            if value in game[i][cellColumn].domain:
                cellRowDomain.remove(value)
                #print celldomain
    for j in range(0,N):
        if((game[cellRow][j].value ==-1) and game[cellRow][j].domain !=None):
            cellColDomain = game[cellRow][j].domain
            if value in game[cellRow][j].domain:
                cellColDomain.remove(value)
        
    x1 = (cellRow/M) * M
    y1 = (cellColumn/K) * K
    for i in range(x1, x1+M):
        for j in range(y1, y1+K):
            if(i >= 0 and i < N and j >= 0 and j < N):
                #if game[i][j].id != cell.id:
                    #if(cellRow == 0 and cellColumn==7 ):
                        #print ("x1 is: ",x1," y1 is: ",y1)
                        #print("i is: ",i," j is: ",j)
                if((game[i][j].value ==-1)  and game[i][j].domain!=None):
                    cellGridDomain = game[i][j].domain
                    if value in game[i][j].domain:
                        cellGridDomain.remove(value)
   
def addDomain(cell,game,value):
    global M,K,N
    #celldomain = cell.domain
    cellRow = cell.row
    cellColumn = cell.column
    #if(cell.value==-1 and celldomain!=None):
        #print celldomain
    for i in range(0,N):
        if((game[i][cellColumn].value ==-1)  and game[i][cellColumn].domain !=None ):
            cellRowDomain = game[i][cellColumn].domain
            if value not in game[i][cellColumn].domain:
                cellRowDomain.append(value)
                #print celldomain
    for j in range(0,N):
        if((game[cellRow][j].value ==-1) and game[cellRow][j].domain !=None):
            cellColDomain = game[cellRow][j].domain
            if value not in game[cellRow][j].domain:
                cellColDomain.append(value)
        
    x1 = (cellRow/M) * M
    y1 = (cellColumn/K) * K
    for i in range(x1, x1+M):
        for j in range(y1, y1+K):
            if(i >= 0 and i < N and j >= 0 and j < N):
                #if game[i][j].id != cell.id:
                    #if(cellRow == 0 and cellColumn==7 ):
                        #print ("x1 is: ",x1," y1 is: ",y1)
                        #print("i is: ",i," j is: ",j)
                if((game[i][j].value ==-1) and game[i][j].domain !=None):
                    cellGridDomain = game[i][j].domain
                    if value not in game[i][j].domain:
                        cellGridDomain.append(value)
    
    
def getminConflictValue(game,cell):
    cellRow = cell.row
    cellColumn = cell.column
    minValueConflicts = float('inf')
    requiredValueForCell = None
    for i in cell.domain:
        cell.value = i
        conflicts = findconflicts(game,cell)
        cell.value = -1
        if conflicts <= minValueConflicts:
            requiredValueForCell = i
    return requiredValueForCell
        
        
def findconflicts(game,cell):
    global M,K,N,conflictsDictionary
    #celldomain = cell.domain
    cellRow = cell.row
    cellColumn = cell.column
    value = cell.value
    conflictCountLocal = 0
    
    #if(cell.value==-1 and celldomain!=None):
        #print celldomain
    for i in range(0,N):
        #print("conflict calculations: ","row is: ",cellRow," column is: ",cellColumn)
        if((game[i][cellColumn].value ==-1)  and game[i][cellColumn].domain !=None ):
            cellRowDomain = game[i][cellColumn].domain
            if value in game[i][cellColumn].domain:
                conflictCountLocal += 1
                
                #print celldomain
    for j in range(0,N):
        if((game[cellRow][j].value ==-1) and game[cellRow][j].domain !=None):
            cellColDomain = game[cellRow][j].domain
            if value in game[cellRow][j].domain:
                conflictCountLocal += 1
        
    x1 = (cellRow/M) * M
    y1 = (cellColumn/K) * K
    for i in range(x1, x1+M):
        for j in range(y1, y1+K):
            if(i >= 0 and i < N and j >= 0 and j < N):
                #if game[i][j].id != cell.id:
                    #if(cellRow == 0 and cellColumn==7 ):
                        #print ("x1 is: ",x1," y1 is: ",y1)
                        #print("i is: ",i," j is: ",j)
                if((game[i][j].value ==-1)  and game[i][j].domain!=None):
                    cellGridDomain = game[i][j].domain
                    if value in game[i][j].domain:
                        conflictCountLocal += 1
    
    return conflictCountLocal 

def prepPriorityQueue(game):
    global sudokuQueue 
    sudokuQueue = PriorityQueue()
    for i in range(0,N):
        for j in range(0,N):
            if((not game[i][j].domain==None) and (len(game[i][j].domain )!=0)):
                sudokuQueue.put((len(game[i][j].domain),game[i][j]))
def readInput(filename):
    global N, M, K, board
    domain = []
    id = 0
    if len(board) != 0:
        return
    board = [[Cell for j in range(N)] for i in range(N)]
    file = open(filename, 'r')
    i = row = 0
    for s in file:
        cell = []
        if i == 0:
            boardConfig = s.split(',')
            N = int(boardConfig[0])
            M = int(boardConfig[1])
            K = int(boardConfig[2].split(';')[0])
            for x in range(1, N + 1):
                domain.append(x)
            i += 1
        else:
            line = s.split(';')[0].split(',')
            column = 0
            for value in line:
                id += 1
                if (value == '-'):
                    c = Cell(id, -1, row, column, False, copy.deepcopy(domain))
                else:
                    c = Cell(id, int(value), row, column, True, None)
                cell.append(c)
                column += 1
            board.append(cell)
            row += 1
    file.close()
    return copy.deepcopy(board)

class Cell:
    def __init__(self, id, value, row, column, readOnly, domain):
        self.id = id
        self.value = value
        self.row = row
        self.column = column
        self.readOnly = readOnly
        self.domain = domain

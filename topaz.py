# In this code, the notation
#   Operator 3: CCNOT(3[1], 1[1], 0)
# will translate to
#   Operator 3 - third operator in the circuit
#   CCNOT(3[1], 1[1], 0) - Apply operator CCNOT.
#   3[1], 1[1], 0 - 3rd and 1st input bits as triggers and 0th input bit as target
#   3[1] - trigger state for 3rd input bit is state 1
#   1[1] - trigger state for 1st input bit is state 1

from numpy import kron, array, identity, ones, zeros, dot, matrix
from math import sqrt, ceil, log, pi, floor
import matplotlib.pyplot as plotter
import logging

print
"hello world!"

maxPatternLength = 2
oneBitAdderBitCount = int(ceil(log(maxPatternLength, 2))) + 1

# translates to maximum pattern length of (2^registerBitCount) characters
registerBitCount = 2

# number of bits alloted for representing a pattern
patternBitCount = 1  # translates to patterns of length 2 at most

# define logging tool
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


# ========= CLASSES ============================================================

class Symbol(object):
    """A class for representing an alphabet, numerical or special character."""

    # attributes:
    #   character - the alphabet, numerical or special character represented by the symbol
    def __init__(self, character):
        super(Symbol, self).__init__()
        self.character = character

    def getCharacter(self):
        return self.character

    # returns the assigned binary code for a symbol
    def getBinaryCode(self):
        if self.character == 'a':
            return [0, 0]
        elif self.character == 'c':
            return [1, 0]
        elif self.character == 't':
            return [0, 1]
        elif self.character == 'g':
            return [1, 1]


class Alphabet(object):
    """A class for representing a finite alphabet."""

    # attributes:
    def __init__(self):
        super(Alphabet, self).__init__()

        # alphabet symbols
        self.symbols = [Symbol('a'), Symbol('c'), Symbol('t'), Symbol('g')]

        # corresponding binary code for each alphabet symbol
        self.codes = [[0, 0], [1, 0], [0, 1], [1, 1]]

        # the number of bits used in binary number for representing a symbol
        self.codeBitCount = 2

        # the number of alphabet symbols
        self.size = len(self.symbols)

    def getSymbols(self):
        return self.symbols

    def getCodeBitCount(self):
        return self.codeBitCount

    def getBinaryCode(self, character):
        for i in xrange(0, self.size):
            if self.symbols[i].character == character:
                return self.codes[i]

    def toString(self):
        text = 'Sigma = {'
        for symbol in self.symbols:
            text = text + symbol + ', '
        text = text + '}'
        return text


class Bit(object):
    """A class for representing a qubit. Its state is either a 0 or a 1."""

    # attributes:
    #   state - either the state 0 or 1
    def __init__(self, state):
        super(Bit, self).__init__()
        self.state = state

    # This method returns the state of this bit which is either a 0 or a 1.
    def getState(self):
        return self.state

    # This is a method which returns the basis vector representation of
    # the state of this bit. Refer to ket0() and ket1() method.
    def getBasisVector(self):
        if self.state == 0:
            return ket0()
        else:
            return ket1()

    # This is a method for adding this bit with an input bit and an input carry bit.
    def add(self, carryBit, bit):
        # print 'BEFORE: Bit.add: carryBit.value =',carryBit.state,', self.value =',self.state,', bit.value =',bit.state
        additionResult = FullAdder().add(carryBit, self, bit)

        # print 'AFTER: Bit.add: carryBit.value =',additionResult.carryBitState,', self.value =',additionResult.sumBitState,', bit.value =',bit.state
        self.state = additionResult.sumBit.getState()

        # return a bit with state equal to carry value resulting from addition
        return Bit(additionResult.carryBit.getState())

    # This is a method which checks if the value(state) of this bit is greater
    # than the value(state) of an input bit. It returns the value True if
    # self > otherBit and False otherwise.
    def isGreaterThan(self, otherBit):
        # print 'Bit.isGreaterThan: self.state:',self.state, ' otherBit.state',otherBit.state
        result = BitComparator().compare(self, otherBit)
        # print 'Bit.isGreaterThan: result:',result
        if result == '11':
            return True
        else:
            return False

    # This method returns the complement bit of this bit. Bit(0) -> Bit(1), Bit(1) -> Bit(0)
    def getComplement(self):
        if self.state == 0:
            return Bit(1)
        else:
            return Bit(0)

    def toString(self):
        return str(self.state)


class BinaryNumber(object):
    """This is a class for representing a binary number. A BinaryNumber is composed of a list of Bit objects."""

    # attributes:
    #   bitValues - a list of bit values representing the value of the binary number
    def __init__(self, bitValues):
        super(BinaryNumber, self).__init__()

        # Note that the number of bits for a binary number is limited to the value of the global variable registerBitCount.
        self.bits = zeroBits()
        for (index, bitState) in enumerate(bitValues):
            self.bits[index].state = bitState

    # This method returns the state of the bits composing this binary number.
    # Output:
    #   state - a list of 0s and 1s corresponding to this BinaryNumber
    def getState(self):
        state = []
        for bit in self.bits:
            state.append(bit.getState())
        return state

    # This method returns the decimal number representation of this binary number.
    def getDecimalNumberFormat(self):
        decimalNumber = 0
        for index, bit in enumerate(self.bits):
            if bit.getState() == 1:
                decimalNumber = decimalNumber + pow(2, index)
        return decimalNumber

    # This method   adds to this binary number an input binary number. The overflow carry value gets lost in the process.
    # TODO: add extra bit for overflow value
    def add(self, otherBinaryNumber):
        # print 'Add:',self.toString(),' +',otherBinaryNumber.toString()
        carryBit = Bit(0)
        for selfBit, otherBit in zip(self.bits, otherBinaryNumber.bits):
            # Note that the resulting sum bit value is put into the current Bit in self.bits and the resulting carry bit value is returned as result of add().
            carryBit = selfBit.add(carryBit, otherBit)

    # This method subtracts from this binary number an input binary number.
    def subtract(self, otherBinaryNumber):
        # Convert the input binary number into 2's complement representation to get its negative value.
        otherTwosComplement = otherBinaryNumber.getTwosComplement()
        # Add to this binary number the negative value of the input binary number. The result of the operation is saved into self.bits.
        self.add(otherTwosComplement)

    # This method returns the complement of this binary number as another binary number.
    def getComplement(self):
        complement = []
        for bit in self.bits:
            complement.append(bit.getComplement().getState())
        return BinaryNumber(complement)

    # This method returns this binary number in 2's-complement notation.
    def getTwosComplement(self):
        complementBinaryNumber = self.getComplement()
        bits = zeros(registerBitCount, dtype=int)
        bits[0] = 1
        oneBinaryNumber = BinaryNumber(bits)
        complementBinaryNumber.add(oneBinaryNumber)
        return complementBinaryNumber

    # This method determines if this binary number is greater than an input binary number.
    def isGreaterThan(self, otherBinaryNumber):
        reversedSelfBitList = self.bits
        reversedSelfBitList.reverse()
        reversedOtherBitList = otherBinaryNumber.bits
        reversedOtherBitList.reverse()

        for (bit, otherBit) in zip(reversedSelfBitList, reversedOtherBitList):
            if bit.isGreaterThan(otherBit):
                return True
            elif otherBit.isGreaterThan(bit):
                return False
        return False

    def toString(self):
        text = '['
        for bit in self.bits:
            text = text + bit.toString()
        text = text + ']'
        return text


class DecimalNumber(object):
    """This class abstracts a decimal number. It provides utility method for converting to binary notation."""

    def __init__(self, value):
        self.value = value

    # This method returns the value of this decimal number.
    def getValue(self):
        return self.value

    # This method returns this decimal number's binary number format.
    def getBinaryNumberFormat(self):
        decimalNumber = int(self.value)
        bitList = []
        for x in xrange(0, registerBitCount):
            i = registerBitCount - x - 1
            if decimalNumber >= pow(2, i):
                bitList.append(1)
                decimalNumber = decimalNumber - pow(2, i)
            else:
                bitList.append(0)
        bitList.reverse()
        return BinaryNumber(bitList)


class BitComparator(object):
    """This class represents a bit-to-bit binary comparator."""

    # attributes:
    #   matrix - this represents the effect of the operator on two input bits to compare
    def __init__(self):
        super(BitComparator, self).__init__()

        # the first 2 multiplicands correspond to bits to compare; last two correspond to auxiliary bits
        matrixDimension = pow(2, 1) * pow(2, 1) * pow(2, 1) * pow(2, 1)
        self.matrix = identity(matrixDimension, dtype=int)

        # Operator 1: CNOT(3[1],1)
        dummyOperator = identity(matrixDimension, dtype=int)
        for bit2 in xrange(0, 1):
            for bit0 in xrange(0, 1):
                decimalValue0 = pow(2, 3) * 1 + pow(2, 2) * bit2 + pow(2, 1) * 0 + pow(2, 0) * bit0
                decimalValue1 = pow(2, 3) * 1 + pow(2, 2) * bit2 + pow(2, 1) * 1 + pow(2, 0) * bit0
                dummyOperator[decimalValue0][decimalValue1] = 1
                dummyOperator[decimalValue1][decimalValue0] = 1
                dummyOperator[decimalValue0][decimalValue0] = 0
                dummyOperator[decimalValue1][decimalValue1] = 0
        self.matrix = dot(self.matrix, dummyOperator)

        # Operator 2: CNOT(2[1],1)
        dummyOperator = identity(matrixDimension, dtype=int)
        for bit3 in xrange(0, 1):
            for bit0 in xrange(0, 1):
                decimalValue0 = pow(2, 3) * bit3 + pow(2, 2) * 1 + pow(2, 1) * 0 + pow(2, 0) * bit0
                decimalValue1 = pow(2, 3) * bit3 + pow(2, 2) * 1 + pow(2, 1) * 1 + pow(2, 0) * bit0
                dummyOperator[decimalValue0][decimalValue1] = 1
                dummyOperator[decimalValue1][decimalValue0] = 1
                dummyOperator[decimalValue0][decimalValue0] = 0
                dummyOperator[decimalValue1][decimalValue1] = 0
        self.matrix = dot(dummyOperator, self.matrix)

        # Operator 3: CCNOT(3[1], 1[1], 0)
        dummyOperator = identity(matrixDimension, dtype=int)
        for bit2 in xrange(0, 1):
            decimalValue0 = pow(2, 3) * 1 + pow(2, 2) * bit2 + pow(2, 1) * 1 + pow(2, 0) * 0
            decimalValue1 = pow(2, 3) * 1 + pow(2, 2) * bit2 + pow(2, 1) * 1 + pow(2, 0) * 1
            dummyOperator[decimalValue0][decimalValue1] = 1
            dummyOperator[decimalValue1][decimalValue0] = 1
            dummyOperator[decimalValue0][decimalValue0] = 0
            dummyOperator[decimalValue1][decimalValue1] = 0
        self.matrix = dot(dummyOperator, self.matrix)

    # This method compares to input bits bit1 and bit2. It returns a 2-bit code as result of the comparison.
    # The first bit serves as an equivalence indicator and the second bit serves as a relation indicator.
    # The following codes are the possible values of the bit comparator and their meanings:
    # 00 -> bit1 == bit2, 10 -> bit1 < bit2, 11 -> bit1 > bit2.
    # TODO: Write in Tex the circuit for 1-bit BitComparator operator
    def compare(self, bit1, bit2):
        basisVector = kron(bit1.getBasisVector(), kron(bit2.getBasisVector(), kron(ket0(), ket0())))
        basisVector = dot(self.matrix, basisVector).tolist()
        indexOf1 = basisVector.index(1)

        # identify value for first auxiliary bit
        equivalenceIndicator = '0'
        startingIndices = [(index * 4) for index in range(4)]
        if indexOf1 in startingIndices or (indexOf1 - 1) in startingIndices:
            equivalenceIndicator = '0'
        elif (indexOf1 - 2) in startingIndices or (indexOf1 - 3) in startingIndices:
            equivalenceIndicator = '1'

        # identify value for second auxiliary bit
        relationIndicator = '0'
        if indexOf1 % 2 == 1:
            relationIndicator = '1'

        return equivalenceIndicator + relationIndicator


# This class represents a full adder for 2 input bits, 1 carry bit, 2 ancilla bit.
class FullAdder(object):
    """This is a class for representing a quantum full adder. It makes use of
        a one-bit adder in computing for the total sum. Refer to OneBitAdder
        class."""

    # attributes:
    #   operator - the matrix representation of a full adder for two bits
    def __init__(self):
        super(FullAdder, self).__init__()

        # TODO: Put into TeX the circuit diagram for the quantum full adder.
        matrixDimension = pow(2, 1) * pow(2, 1) * pow(2, 1) * pow(2, 1) * pow(2, 1)
        self.matrix = identity(matrixDimension, dtype=int)

        # Operator 1: CCNOT(3[1],2[1],1)
        dummyOperator = identity(matrixDimension, dtype=int)
        for bit4 in xrange(0, 2):
            for bit0 in xrange(0, 2):
                decimal_value_0 = pow(2, 4) * bit4 + pow(2, 3) * 1 + pow(2, 2) * 1 + pow(2, 1) * 0 + pow(2, 0) * bit0
                decimal_value_1 = pow(2, 4) * bit4 + pow(2, 3) * 1 + pow(2, 2) * 1 + pow(2, 1) * 1 + pow(2, 0) * bit0
                dummyOperator[decimal_value_0][decimal_value_1] = 1
                dummyOperator[decimal_value_1][decimal_value_0] = 1
                dummyOperator[decimal_value_0][decimal_value_0] = 0
                dummyOperator[decimal_value_1][decimal_value_1] = 0
        self.matrix = dot(dummyOperator, self.matrix)

        # Operator 2: CNOT(3[1],2)
        dummyOperator = identity(matrixDimension, dtype=int)
        for bit4 in xrange(0, 2):
            for bit1 in xrange(0, 2):
                for bit0 in xrange(0, 2):
                    decimal_value_0 = pow(2, 4) * bit4 + pow(2, 3) * 1 + pow(2, 2) * 0 + pow(2, 1) * bit1 + pow(2,
                                                                                                                0) * bit0
                    decimal_value_1 = pow(2, 4) * bit4 + pow(2, 3) * 1 + pow(2, 2) * 1 + pow(2, 1) * bit1 + pow(2,
                                                                                                                0) * bit0
                    dummyOperator[decimal_value_0][decimal_value_1] = 1
                    dummyOperator[decimal_value_1][decimal_value_0] = 1
                    dummyOperator[decimal_value_0][decimal_value_0] = 0
                    dummyOperator[decimal_value_1][decimal_value_1] = 0
        self.matrix = dot(dummyOperator, self.matrix)

        # Operator 3: CCNOT(4[1],2[1],0)
        dummyOperator = identity(matrixDimension, dtype=int)
        for bit3 in xrange(0, 2):
            for bit1 in xrange(0, 2):
                decimal_value_0 = pow(2, 4) * 1 + pow(2, 3) * bit3 + pow(2, 2) * 1 + pow(2, 1) * bit1 + pow(2, 0) * 0
                decimal_value_1 = pow(2, 4) * 1 + pow(2, 3) * bit3 + pow(2, 2) * 1 + pow(2, 1) * bit1 + pow(2, 0) * 1
                dummyOperator[decimal_value_0][decimal_value_1] = 1
                dummyOperator[decimal_value_1][decimal_value_0] = 1
                dummyOperator[decimal_value_0][decimal_value_0] = 0
                dummyOperator[decimal_value_1][decimal_value_1] = 0
        self.matrix = dot(dummyOperator, self.matrix)

        # Operator 4: CNOT(4[1],2)
        dummyOperator = identity(matrixDimension, dtype=int)
        for bit3 in xrange(0, 2):
            for bit1 in xrange(0, 2):
                for bit0 in xrange(0, 2):
                    decimal_value_0 = pow(2, 4) * 1 + pow(2, 3) * bit3 + pow(2, 2) * 0 + pow(2, 1) * bit1 + pow(2,
                                                                                                                0) * bit0
                    decimal_value_1 = pow(2, 4) * 1 + pow(2, 3) * bit3 + pow(2, 2) * 1 + pow(2, 1) * bit1 + pow(2,
                                                                                                                0) * bit0
                    dummyOperator[decimal_value_0][decimal_value_1] = 1
                    dummyOperator[decimal_value_1][decimal_value_0] = 1
                    dummyOperator[decimal_value_0][decimal_value_0] = 0
                    dummyOperator[decimal_value_1][decimal_value_1] = 0
        self.matrix = dot(dummyOperator, self.matrix)

        # Operator 5: CNOT(1[1],0)
        dummyOperator = identity(matrixDimension, dtype=int)
        for bit4 in xrange(0, 2):
            for bit3 in xrange(0, 2):
                for bit2 in xrange(0, 2):
                    decimal_value_0 = pow(2, 4) * bit4 + pow(2, 3) * bit3 + pow(2, 2) * bit2 + pow(2, 1) * 1 + pow(2,
                                                                                                                   0) * 0
                    decimal_value_1 = pow(2, 4) * bit4 + pow(2, 3) * bit3 + pow(2, 2) * bit2 + pow(2, 1) * 1 + pow(2,
                                                                                                                   0) * 1
                    dummyOperator[decimal_value_0][decimal_value_1] = 1
                    dummyOperator[decimal_value_1][decimal_value_0] = 1
                    dummyOperator[decimal_value_0][decimal_value_0] = 0
                    dummyOperator[decimal_value_1][decimal_value_1] = 0
        self.matrix = dot(dummyOperator, self.matrix)

    # This is a method for adding two bit values. e.g. 1 + 1 = 0 carry 1
    # Input:
    #   1. carryBitValue - the state of a carry bit
    #   2. bitValue1 - the state of the first bit
    #   3. bitValue2 - the state of the second bit.
    # Output:
    #   1. additionResult - an AdditionResult object
    def add(self, carryBit, bit1, bit2):

        # prepare the basis vector representation of the state of the register
        registerBasisVector = kron(carryBit.getBasisVector(),
                                   kron(bit1.getBasisVector(),
                                        kron(bit2.getBasisVector(),
                                             kron(ket0(), ket0()))))

        # perform binary addition by performing dot multiplication on the
        # basis vector representing the register's state and the full adder's
        # matrix representation
        registerBasisVector = dot(self.matrix, registerBasisVector).tolist()
        indexOf1 = registerBasisVector.index(1)

        # store the resulting bit from the bit-to-bit addition into sumBitValue
        sumBitValue = 0
        for startIndex in [0, 8, 16, 24]:
            if indexOf1 in range(startIndex, startIndex + 8):
                if indexOf1 in range(startIndex, startIndex + 4):
                    sumBitValue = 0
                else:
                    sumBitValue = 1

        # store the resulting carry bit from the bit-to-bit addition into
        # carryBitValue
        carryBitValue = 0
        if indexOf1 % 2 == 1:
            carryBitValue = 1

        # store the resulting sum bit and carry bit into an AdditionResult
        # object
        additionResult = AdditionResult(Bit(sumBitValue), Bit(carryBitValue))
        return additionResult


class AdditionResult(object):
    """This is a class for representing the result of a binary addition using
        class FullAdder."""

    # attributes:
    #   sumBit - the resulting sum bit from binary addition of class type Bit
    #   carryBit - the resulting carry bit from binary addition of class type Bit
    def __init__(self, sumBit, carryBit):
        super(AdditionResult, self).__init__()
        self.sumBit = sumBit
        self.carryBit = carryBit

    def getSumBit(self):
        return self.sumBit

    def getCarryBit(self):
        return self.carryBit


class U_Loc(object):
    """This is an operator for identifying the index of first occurrence of each symbol in text T in pattern P."""

    # attributes:
    #   text - the input text
    #   pattern - the input pattern
    #   P_Sym - the set of unique symbols in P
    #   P_Loc - the set of indices in P of symbols in P_Sym
    def __init__(self, text, pattern, P_Sym, P_Loc):
        super(U_Loc, self).__init__()
        textLength = len(text)
        patternLength = len(pattern)
        matrixDimension = textLength * textLength
        self.matrix = zeros((matrixDimension, matrixDimension), dtype=int)
        for textIndex in xrange(0, textLength):
            if text[textIndex] in P_Sym:
                P_Sym_index = P_Sym.index(text[textIndex])  # the index of the symbol T[i] in set P_Sym
                P_Loc_index = P_Loc[P_Sym_index]  # the index of the symbol T[i] in pattern P
                self.matrix[textIndex * textLength][textIndex * textLength + P_Loc_index] = 1
                self.matrix[textIndex * textLength + P_Loc_index][textIndex * textLength] = 1
                for i in xrange(0, patternLength):
                    if i != 0 and i != P_Loc_index:
                        self.matrix[textIndex * textLength + i][textIndex * textLength + i] = 1
                for i in xrange(patternLength, pow(2, registerBitCount)):
                    self.matrix[(textIndex * textLength) + i][(textIndex * textLength) + i] = 1
            else:
                for i in xrange(0, pow(2, registerBitCount)):
                    self.matrix[textIndex * textLength + i][textIndex * textLength + i] = 1
        print
        'self.matrix:'
        print
        self.matrix

    # This method applies the U_Loc operator to an input vector of appropriate dimensions
    def apply(self, phase1State):
        if len(phase1State.getBasisVector()) == self.matrix.shape[0] and len(phase1State.getBasisVector()) ==
                self.matrix.shape[1]:
            resultingBasisVector = dot(self.matrix, phase1State.getBasisVector())
            patternBasisVector = resultingBasisVector[phase1State.getTextIndexStateDecimal() * pow(2,
                                                                                                   patternBitCount): phase1State.getTextIndexStateDecimal() * pow(
                2, patternBitCount) + pow(2, patternBitCount)]
            phase1State.updatePatternIndex(dot(self.matrix, phase1State.getBasisVector()))
        else:
            print
            'The input vector\'s dimension does not match the dimension of the U_Loc operator\'s matrix.'


class Phase1State(object):
    """This class represents the quantum state of register for representing text index and another register for representing pattern index."""

    # Initialize Phase1State object.
    # Input:
    #   amplitude - a floating point number representing total amplitude of state representing text index and state representing pattern index
    #   textIndexState - a BinaryNumber corresponding to state representing a text index (0 -> N-1)
    #   patternIndexState - a BinaryNumber corresponding to state representing a pattern index (0 -> M-1)
    def __init__(self, amplitude, textIndexState, patternIndexState):
        super(Phase1State, self).__init__()
        self.amplitude = amplitude
        self.textIndexState = textIndexState
        self.patternIndexState = patternIndexState

    def getAmplitude(self):
        return self.amplitude

    def getTextIndexState(self):
        return self.textIndexState

    # This method returns the decimal number representation of the binary number corresponding to the state of the text index register.
    def getTextIndexStateDecimal(self):
        return self.textIndexState.getDecimalNumberFormat()

    # This method returns the basis vector representation of the binary number corresponding to the state of the text index register.
    def getTextIndexStateBasisVector(self):
        vector = zeros(pow(2, registerBitCount), dtype=int).tolist()
        vector[self.getTextIndexStateDecimal()] = 1
        return vector

    def getPatternIndexState(self):
        return self.patternIndexState

    # This method returns the decimal number representation of the binary number corresponding to the state of the pattern index register.
    def getPatternIndexStateDecimal(self):
        return self.patternIndexState.getDecimalNumberFormat()

    # This method returns the basis vector representation of the binary number corresponding to the state of the pattern index register.
    def getPatternIndexStateBasisVector(self):
        vector = zeros(pow(2, registerBitCount), dtype=int).tolist()
        #        vector = zeros(pow(2,patternBitCount),dtype=int).tolist()
        vector[self.getPatternIndexStateDecimal()] = 1
        return vector

    # Update the pattern index associated with this Phase1State. This method is called when applying operator U_Loc to an initial Phase1State to reflect the resulting Phase1State.
    def updatePatternIndex(self, newBasisVector):
        patternBasisVector = newBasisVector[self.getTextIndexStateDecimal() * pow(2,
                                                                                  registerBitCount): self.getTextIndexStateDecimal() * pow(
            2, registerBitCount) + pow(2, registerBitCount)]
        self.patternIndexState = DecimalNumber(list(patternBasisVector).index(1)).getBinaryNumberFormat()

    # This method computes i - Phi(T[i]) where Phi(T[i]) is the index of first occurrence in pattern P of symbol T[i].
    # Output:
    #   self.patternIndexState = i - Phi(T[i]), replaces the state self.patternIndexState with the value i - Phi(T[i]) in BinaryNumber format
    def computeLocation(self):
        textIndexBinaryNumber = BinaryNumber(self.textIndexState.getState())
        textIndexBinaryNumber.subtract(self.patternIndexState)
        self.patternIndexState = textIndexBinaryNumber

    # This method returns the kron of the basis vector for the text index state and the basis vector for the pattern index state.
    def getBasisVector(self):
        return kron(self.getTextIndexStateBasisVector(), self.getPatternIndexStateBasisVector())

    def toString(self):
        return 'amplitude:' + str(self.amplitude) + ', text index state:|' + str(
            self.getTextIndexStateDecimal()) + ',' + self.textIndexState.toString() + ',' + str(
            self.getTextIndexStateBasisVector()) + ',>, pattern index state:|' + str(
            self.getPatternIndexStateDecimal()) + ',' + self.patternIndexState.toString() + ',' + str(
            self.getPatternIndexStateBasisVector()) + '>, basis vector:' + str(self.getBasisVector())


# Given a superposition state for phase 1 of the algorithm, this method creates a dictionary (hashmap) where the key values are the states |textIndexState - patternIndexState> and the values are the states |textIndexState> associated with their corresponding |textIndexState - patternIndexState> states. |textIndexState - patternIndexState> states with more associated |textIndexState> states will have a higher amplitude and signify a highly probable starting index of pattern P in text T.
def expressSuperpositionInCompressedRepresentation(phase1SuperpositionState, textLength):
    dict = {key: [] for key in xrange(0, textLength * textLength)}
    for phase1State in phase1SuperpositionState:
        dict[phase1State.getPatternIndexStateDecimal()].append(
            (phase1State.getAmplitude(), phase1State.getTextIndexStateDecimal()))

    return dict


# Given a dictionary output from method expressSuperpositionInCompressedRepresentation compute the total probability associated with each key |textIndexState - patternIndexState>. Return a new dictionary with key |textIndexState - patternIndexState> and value total probability.
def computeTotalProbabilityForEachState(stateToAmplitudeMap):
    dict = {}
    for key in stateToAmplitudeMap.keys():
        probability = 0.0
        for amplitudeStateTuple in stateToAmplitudeMap[key]:
            probability = probability + (amplitudeStateTuple[0] * amplitudeStateTuple[0])
        dict[key] = probability
    return dict


# Given a map in which keys correspond to state index and values correspond to probabilities, plot a probability bar chart for all states.
def displayProbabilityChart(dict):
    #    dict = {'hello': 3, 'how': 4, 'yes': 10, 'you': 11, 'days': 10, 'are': 20, 'ago': 11}

    plotter.title('Probability of occurence for each index in T')
    plotter.ylabel('Probability')
    plotter.xlabel('Index')
    plotter.xticks(range(len(dict)), list(dict.keys()))
    width = 0.05
    plotter.bar(range(len(dict)), dict.values(), width, color="black", align="center")
    plotter.grid(True)
    plotter.ylim(0.0, 1.0)
    plotter.xlim(0, 16)
    plotter.show()


# ================================== Tests ======================================

# ================= Symbol ===================
print
'>> Testing class Symbol.'
codes = [[0, 0], [1, 0], [0, 1], [1, 1]]
for index, character in enumerate(['a', 'c', 't', 'g']):
    print
    'Character:', character
    print
    'Code:', codes[index]
    symbol = Symbol(character)
    if (symbol.getCharacter() == character):
        print
        'Symbol character is correct.'
    else:
        print
        'Symbol character is wrong.'
    if (symbol.getBinaryCode() == codes[index]):
        print
        'Symbol code is correct.'
    else:
        print
        'Symbol code is wrong'

# ================= Alphabet ===================
print
''
print
'>> Testing class Alphabet.'
symbols = [Symbol('a'), Symbol('c'), Symbol('t'), Symbol('g')]
alphabet = Alphabet()
if (alphabet.codes == codes):
    print
    'Alphabet codes are correct.'
else:
    print
    'Alphabet codes are wrong.'
for index, symbol in enumerate(symbols):
    if (alphabet.getBinaryCode(symbol.getCharacter()) == codes[index]):
        print
        'Binary code for character', symbol.getCharacter(), ' is correct.'
    else:
        print
        'Binary code for character', symbol.getCharacter(), ' is wrong.'


def ket0():
    return array([1, 0])


def ket1():
    return array([0, 1])


# This method returns a list of Bit objects in state 0. The list takes the
# global variable registerBitCount for its size.
def zeroBits():
    bits = []
    for i in xrange(0, registerBitCount):
        bits.append(Bit(0))
    return bits


# ================= Bit ===================
print
''
print
'>>>> Testing class Bit.'
print
'>> Testing getter methods.'
bit = Bit(0)
if (bit.getState() == 0):
    print
    'Bit 0 state is correct.'
else:
    print
    'Bit 0 state is incorrect.'
if (bit.getBasisVector().all() == ket0().all()):
    print
    'Bit 0 basis vector is correct.'
else:
    print
    'Bit 0 basis vector is incorrect.'

bit = Bit(1)
if (bit.getState() == 1):
    print
    'Bit 1 state is correct.'
else:
    print
    'Bit 1 state is incorrect.'
if (bit.getBasisVector().all() == ket1().all()):
    print
    'Bit 1 basis vector is correct.'
else:
    print
    'Bit 1 basis vector is incorrect.'

print
'>> Testing Bit.add()'
for bitState in [0, 1]:
    for dummyBitState in [0, 1]:
        bit = Bit(bitState)
        dummyBit = Bit(dummyBitState)
        carryBit = Bit(0)
        print
        bit.getState(), '+', dummyBit.getState(), '=',
        carryBit = bit.add(carryBit, dummyBit)
        print
        bit.getState(), 'carry', carryBit.getState()

print
'>> Testing Bit.isGreaterThan()'
for bitState in [0, 1]:
    for dummyBitState in [0, 1]:
        bit = Bit(bitState)
        dummyBit = Bit(dummyBitState)
        print
        bit.getState(), '>', dummyBit.getState(), ':', (bit.isGreaterThan(dummyBit) == True)

print
''
print
'>> Testing Bit.getComplement()'
for bitState in [0, 1]:
    bit = Bit(bitState)
    print
    bit.getState(), 'complement is', bit.getComplement().getState()

# ================= BitComparator ==================
print
''
print
'>>>> Testing class BitComparator.'
print
'>> Testing BitComparator.compare()'
for bit0State in [0, 1]:
    for bit1State in [0, 1]:
        comparator = BitComparator()
        result = comparator.compare(Bit(bit0State), Bit(bit1State))
        relation = ''
        if result == '00':
            relation = '=='
        elif result == '10':
            relation = '<'
        elif result == '11':
            relation = '>'
        print
        bit0State, relation, bit1State

# ================= BinaryNumber ===================
print
''
print
'>>>> Testing class BinaryNumber.'
print
'>> Testing BinaryNumber.getComplement()'
for bit0State in [0, 1]:
    for bit1State in [0, 1]:
        bin = BinaryNumber([bit0State, bit1State])
        binComplement = bin.getComplement()
        print
        bin.toString(), 'complement', binComplement.toString()

print
''
print
'>> Testing BinaryNumber.getTwosComplement().'
for bit0State in [0, 1]:
    for bit1State in [0, 1]:
        bin = BinaryNumber([bit0State, bit1State])
        binTwosComplement = bin.getTwosComplement()
        print
        bin.toString(), '2\'s complement', binTwosComplement.toString()

print
''
print
'>> Testing BinaryNumber.add().'
for binBit1 in [0, 1]:
    for binBit0 in [0, 1]:
        for dummyBinBit1 in [0, 1]:
            for dummyBinBit0 in [0, 1]:
                bin = BinaryNumber([binBit0, binBit1])
                dummyBin = BinaryNumber([dummyBinBit0, dummyBinBit1])
                print
                bin.toString(), '+', dummyBin.toString(), '=',
                bin.add(dummyBin)
                print
                bin.toString()

print
''
print
'>> Testing BinaryNumber.substract().'
for binBit0 in [0, 1]:
    for dummyBinBit0 in [0, 1]:
        bin = BinaryNumber([binBit0, 0])
        dummyBin = BinaryNumber([dummyBinBit0, 0])
        print
        bin.toString(), '-', dummyBin.toString(), '=',
        bin.subtract(dummyBin)
        print
        bin.toString()

print
''
print
'>> Testing BinaryNumber.isGreaterThan().'
for bin1bit0State in [0, 1]:
    for bin1bit1State in [0, 1]:
        for bin2bit0State in [0, 1]:
            for bin2bit1State in [0, 1]:
                binary1 = BinaryNumber([bin1bit0State, bin1bit1State])
                binary2 = BinaryNumber([bin2bit0State, bin2bit1State])
                print
                binary1.getState(), '>', binary2.getState(), ':', (binary1.isGreaterThan(binary2) == True)

print
''
print
'>> Testing BinaryNumber.getDecimalNumberFormat().'
for bin1bit0State in [0, 1]:
    for bin1bit1State in [0, 1]:
        binary1 = BinaryNumber([bin1bit0State, bin1bit1State])
        print
        'Binary number', binary1.getState(), '== Decimal number', binary1.getDecimalNumberFormat()

print
''
print
'>>>> Testing class DecimalNumber.'
print
'>> Testing DecimalNumber.getValue()'
for index in xrange(0, pow(2, registerBitCount)):
    decimalNumber = DecimalNumber(index)
    print
    'Decimal number', index, ':', decimalNumber.getValue()
    print
    'Decimal number', index, 'in binary format:', decimalNumber.getBinaryNumberFormat().toString()

print
''
print
'>>>> Testing class Phase1State.'
phase1SuperpositionState = []
for textIndex in xrange(0, 4):
    phase1SuperpositionState.append(
        Phase1State(1.0, DecimalNumber(textIndex).getBinaryNumberFormat(), DecimalNumber(0).getBinaryNumberFormat()))
print
'>> Testing Phase1State.toString()'

for index in xrange(0, 4):
    # Phase1State.getAmplitude()
    if phase1SuperpositionState[index].getAmplitude() == 1.0:
        print
        'Amplitude for superpositioned state |', phase1SuperpositionState[index].getTextIndexState().getState(), ',', \
        phase1SuperpositionState[index].getPatternIndexState().getState(), '> is', phase1SuperpositionState[
            index].getAmplitude(), 'and is correct.'
    else:
        print
        'Amplitude for superpositioned state |', phase1SuperpositionState[index].getTextIndexState().getState(), ',', \
        phase1SuperpositionState[index].getPatternIndexState().getState(), '> is', phase1SuperpositionState[
            index].getAmplitude(), 'and is wrong.'

    # Phase1State.getTextIndexStateDecimal()
    if phase1SuperpositionState[index].getTextIndexStateDecimal() == index:
        print
        'Decimal number representation', phase1SuperpositionState[
            index].getTextIndexStateDecimal(), 'for text index binary number state', phase1SuperpositionState[
            index].getTextIndexState().getState(), 'is correct.'
    else:
        print
        'Decimal number representation', phase1SuperpositionState[
            index].getTextIndexStateDecimal(), 'for text index binary number state', phase1SuperpositionState[
            index].getTextIndexState().getState(), 'is wrong.'

    # Phase1State.getPatternIndexStateDecimal()
    if phase1SuperpositionState[index].getPatternIndexStateDecimal() == 0:
        print
        'Decimal number representation', phase1SuperpositionState[
            index].getPatternIndexStateDecimal(), 'for pattern index binary number state', phase1SuperpositionState[
            index].getPatternIndexState().getState(), 'is correct.'
    else:
        print
        'Decimal number representation', phase1SuperpositionState[
            index].getPatternIndexStateDecimal(), 'for pattern index binary number state', phase1SuperpositionState[
            index].getPatternIndexState().getState(), 'is wrong.'
    print
    ''

# ========================== Text processing ====================================
# a. Ask for input text. #TODO: Implement IO for input text.
text = 'acgc'
# b. Ask for input pattern. #TODO: Implement IO for input pattern.
pattern = 'gc'
# c. Ask for threshold on number of mismatches. #TODO: Implement IO for input threshold number of mismatches.
threshold = 1  # The maximum number of allowed mismatches between any solution substring of text and the whole pattern
N = len(text)  # The length of the text
M = len(pattern)  # The length of the pattern
P_Sym = list(set(pattern))  # The set of unique symbols in pattern
P_Loc = map((lambda symbol: pattern.index(symbol)),
            P_Sym)  # The location of first occurrence in pattern of symbols in P_Sym
inputParameters = '================== Input parameters >> Text:\'' + text + '\', Pattern:\'' + pattern + '\', Threshold:' + str(
    threshold) + ', P_Sym:' + str(P_Sym) + ', P_Loc:' + str(P_Loc)
logging.info(inputParameters)

# I. Filtering phase ===========================================================
# I.a Preparation of superposition state representing unfiltered search space. =
logging.info('================== Preparing initial superposition state for Phase 1.')
phase1SuperpositionState = []
for textIndex in xrange(0, N):
    phase1SuperpositionState.append(Phase1State(sqrt(1.0 / N), DecimalNumber(textIndex).getBinaryNumberFormat(),
                                                DecimalNumber(0).getBinaryNumberFormat()))
print
''
for textIndex in xrange(0, N):
    print
    phase1SuperpositionState[textIndex].toString()
print
''
logging.info('================== Initial superposition state for Phase 1 prepared.')
logging.info('================== Identifying index of first occurrence in pattern P of each symbol in text T.')
# I.b Identifying location of first occurence of each symbol in text T. ========
locationOperator = U_Loc(text, pattern, P_Sym, P_Loc)
for phase1State in phase1SuperpositionState:
    print
    'State |', phase1State.getTextIndexStateDecimal(), '>|', phase1State.getPatternIndexStateDecimal(), '>:', phase1State.toString()
    locationOperator.apply(phase1State)
    print
    'State |', phase1State.getTextIndexStateDecimal(), '>|', phase1State.getPatternIndexStateDecimal(), '>:', phase1State.toString()
logging.info(
    '================== Identitification of index of first occurrence in pattern P of each symbol in text T is done.')
# I.c Identifying possible starting location of pattern P in text T. ===========
# TODO: Consider negative solutions to textIndex - indexOfFirstOccurrence operation
logging.info('================== Identitification of candidate starting index of pattern P in text T.')
for phase1State in phase1SuperpositionState:
    print
    'State |', phase1State.getTextIndexStateDecimal(), '>|', phase1State.getPatternIndexStateDecimal(), '>: ', phase1State.toString()
    phase1State.computeLocation()
    print
    'State |', phase1State.getTextIndexStateDecimal(), '>|', phase1State.getPatternIndexStateDecimal(), '>: ', phase1State.toString()
logging.info('================== Identitification of candidate starting index of pattern P in text T is done.')

# Code for grouping textIndexStates with same patterIndexState = textIndex - indexOfFirstOccurrence
stateToAmplitudeMap = expressSuperpositionInCompressedRepresentation(phase1SuperpositionState, N)
print
'State to amplitude map:', stateToAmplitudeMap

# Code for computing total probability for each key in stateToAmplitudeMap
stateToProbabilityMap = computeTotalProbabilityForEachState(stateToAmplitudeMap)
print
'State to probability map:', stateToProbabilityMap

# TODO: Write code for displaying probabilities into a bar chart.
displayProbabilityChart(stateToProbabilityMap)

# I.d Measuring state of whole unfiltered search space register to identify
#       possible starting locations of P in T. =================================
# TODO: Write code for measuring state of register during filtering phase.

# II. Verification Phase =======================================================
# II.a Preparation of superposition state representing candidate substrings of T
# TODO: Write code for preparing superposition state for verification phase.
# TODO: Write code for verifying candidate solution indices in text T.
# TODO: Write code for amplifying amplitude of solution states.
# TODO: Write code for measuring state of register during filtering phase.






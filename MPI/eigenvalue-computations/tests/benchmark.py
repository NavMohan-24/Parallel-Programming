from shutil import which
import qutip as qt
import numpy as np

def hamiltonian_ising(N,J,h):
    sx = qt.operators.sigmax()
    sz = qt.operators.sigmaz()
    I = qt.operators.identity(2) 

    H = 0

    for i in range(N-1):
        op_list = [I]*N
        op_list[i] = sz
        op_list[i+1] = sz
        H += -J*qt.tensor(op_list)

    for i in range(N):
        op_list = [I]*N
        op_list[i] = sx
        H += -h* qt.tensor(op_list)
    
    return H


N = 2
J = -1.0
h = 0.1

# Hamiltonian
H = hamiltonian_ising(N,J,h)
print(H)

# collapse operator
kappa = 1
sm = qt.operators.sigmam()
I = qt.operators.identity(2)
c_ops = [qt.tensor([sm if i==j  else I for j in range(N)]) for i in range(N)]
c_ops = [np.sqrt(kappa) * c_op for c_op in c_ops]

# Generating Linbladian
Linbladian = qt.liouvillian(H,c_ops)
# print(Linbladian)

eigenvals, eigenmat = Linbladian.eigenstates(eigvals=16)

for i,vals in enumerate(eigenvals):
    print(f"Eigenvalue {i+1}: ",vals)

for vec in eigenmat:
    print(vec)

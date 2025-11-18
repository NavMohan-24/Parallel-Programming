import matplotlib.pyplot as plt
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


N = 3
J = -1.0
h = 0.1
num_vals = 64


# Hamiltonian
H = hamiltonian_ising(N,J,h)
# print(H)

# collapse operator
kappa = 1
sm = qt.operators.sigmam()
I = qt.operators.identity(2)
c_ops = [qt.tensor([sm if i==j  else I for j in range(N)]) for i in range(N)]
c_ops = [np.sqrt(kappa) * c_op for c_op in c_ops]

# Generating Linbladian
Linbladian = qt.liouvillian(H,c_ops)
# print(Linbladian)

eigenvals, eigenmat = Linbladian.eigenstates(eigvals=num_vals, sort='high')

for i,vals in enumerate(eigenvals):
    print(f"Eigenvalue {i+1}: ",vals)

# for vec in eigenmat:
#     print(vec)

eigenvals.real[np.abs(eigenvals.real) < 1e-12] = 0
eigenvals.imag[np.abs(eigenvals.imag) < 1e-12] = 0
print(eigenvals)

# Load C++ eigenvalues
data10 = np.loadtxt('eigVals10.csv', delimiter=',', skiprows=0)
cpp_eigenvals10 = data10[:, 0] + 1j * data10[:, 1]
cpp_eigenvals10.real[np.abs(cpp_eigenvals10.real) < 1e-12] = 0
cpp_eigenvals10.imag[np.abs(cpp_eigenvals10.imag) < 1e-12] = 0

data20 = np.loadtxt('eigVals20.csv', delimiter=',', skiprows=0)
cpp_eigenvals20 = data20[:, 0] + 1j * data20[:, 1]
cpp_eigenvals20.real[np.abs(cpp_eigenvals20.real) < 1e-12] = 0
cpp_eigenvals20.imag[np.abs(cpp_eigenvals20.imag) < 1e-12] = 0

# data30 = np.loadtxt('eigVals30.csv', delimiter=',', skiprows=0)
# cpp_eigenvals30 = data30[:, 0] + 1j * data30[:, 1]
# cpp_eigenvals30.real[np.abs(cpp_eigenvals30.real) < 1e-12] = 0
# cpp_eigenvals30.imag[np.abs(cpp_eigenvals30.imag) < 1e-12] = 0

data40 = np.loadtxt('eigVals20.csv', delimiter=',', skiprows=0)
cpp_eigenvals40 = data40[:, 0] + 1j * data40[:, 1]
cpp_eigenvals40.real[np.abs(cpp_eigenvals40.real) < 1e-12] = 0
cpp_eigenvals40.imag[np.abs(cpp_eigenvals40.imag) < 1e-12] = 0

# data50 = np.loadtxt('eigVals50.csv', delimiter=',', skiprows=0)
# cpp_eigenvals50 = data50[:, 0] + 1j * data50[:, 1]
# cpp_eigenvals50.real[np.abs(cpp_eigenvals50.real) < 1e-12] = 0
# cpp_eigenvals50.imag[np.abs(cpp_eigenvals50.imag) < 1e-12] = 0

dataf = np.loadtxt('eigValsfull.csv', delimiter=',', skiprows=0)
cpp_eigenvalsf = dataf[:, 0] + 1j * dataf[:, 1]
cpp_eigenvalsf.real[np.abs(cpp_eigenvalsf.real) < 1e-12] = 0
cpp_eigenvalsf.imag[np.abs(cpp_eigenvalsf.imag) < 1e-12] = 0



# print(cpp_eigenvals)


# --- 2. Create the plot ---

fig,ax = plt.subplots(nrows=2, ncols=2,figsize=(10,8))

ax[0,0].scatter(eigenvals.real, eigenvals.imag, color='red', marker='o', label='Exact Spectrum (QuTiP)')
ax[0,0].scatter(cpp_eigenvals10.real, cpp_eigenvals10.imag, color='blue', marker='x', label='Krylov subspace diagonalization')
ax[0,0].title.set_text(r"Krylov subspace dimension $k=10$")
ax[0,0].set_xlim(-4, 1)
ax[0,0].axhline(0, color='black', linewidth=0.5) # Add real axis
ax[0,0].axvline(0, color='black', linewidth=0.5) 
ax[0,0].set_xlabel(r"$\mathit{Re}[\lambda_j]$ (Decay rate)")
ax[0,0].set_ylabel(r"$\mathit{Im}[\lambda_j]$ (Oscillation frequency)")

ax[0,1].scatter(eigenvals.real, eigenvals.imag, color='red', marker='o', label='Exact Spectrum (QuTiP)')
ax[0,1].scatter(cpp_eigenvals20.real, cpp_eigenvals20.imag, color='blue', marker='x', label='Krylov subspace diagonalization')
ax[0,1].title.set_text(r"Krylov subspace dimension $k=20$")
ax[0,1].set_xlim(-4, 1)
ax[0,1].axhline(0, color='black', linewidth=0.5) # Add real axis
ax[0,1].axvline(0, color='black', linewidth=0.5) 
ax[0,1].set_xlabel(r"$\mathit{Re}[\lambda_j]$ (Decay rate)")
ax[0,1].set_ylabel(r"$\mathit{Im}[\lambda_j]$ (Oscillation frequency)")


ax[1,0].scatter(eigenvals.real, eigenvals.imag, color='red', marker='o', label='Exact Spectrum (QuTiP)')
ax[1,0].scatter(cpp_eigenvals40.real, cpp_eigenvals40.imag, color='blue', marker='x', label='Krylov subspace diagonalization')
ax[1,0].title.set_text(r"Krylov subspace dimension $k=40$")
ax[1,0].set_xlim(-4, 1)
ax[1,0].axhline(0, color='black', linewidth=0.5) # Add real axis
ax[1,0].axvline(0, color='black', linewidth=0.5) 
ax[1,0].set_xlabel(r"$\mathit{Re}[\lambda_j]$ (Decay rate)")
ax[1,0].set_ylabel(r"$\mathit{Im}[\lambda_j]$ (Oscillation frequency)")

ax[1,1].scatter(eigenvals.real, eigenvals.imag, color='red', marker='o', label='Exact Spectrum (QuTiP)')
ax[1,1].scatter(cpp_eigenvalsf.real, cpp_eigenvalsf.imag, color='blue', marker='x', label='Krylov subspace diagonalization')
ax[1,1].title.set_text("Full Liouvillian spectrum")
ax[1,1].set_xlim(-4, 1)
ax[1,1].axhline(0, color='black', linewidth=0.5) # Add real axis
ax[1,1].axvline(0, color='black', linewidth=0.5) 
ax[1,1].set_xlabel(r"$\mathit{Re}[\lambda_j]$ (Decay rate)")
ax[1,1].set_ylabel(r"$\mathit{Im}[\lambda_j]$ (Oscillation frequency)")
    
# Add grid lines to all subplots easily using a loop
for a in ax.flat:
    a.grid(True, color='grey', linestyle='--', linewidth=0.5)
    # a.set_facecolor("smoke")
    
    # a.legend()

handles, labels = [], []
for a in ax.flatten():
    for handle, label in zip(*a.get_legend_handles_labels()):
        # Optional: Prevent duplicate labels if a label appears on multiple plots
        if label not in labels:
            handles.append(handle)
            labels.append(label)

plt.suptitle(r"Liouvillian Spectrum of the TFIM ($N=3$, $J=-1.0$, $h=0.1$)", fontsize=14,fontweight='bold')

# --- LAYOUT FIX STARTS HERE ---
# Adjust layout so subplots don’t collide with title or legend
plt.tight_layout(rect=[0, 0.05, 0.95, 0.95])  # leave space at top/right

fig.set_facecolor("whitesmoke")
# Add global legend
fig.legend(
    handles, labels,
    loc='lower center',
    bbox_to_anchor=(0.5, 0.0),  # slightly inset inside figure
    ncol=2,
    fancybox=True,
    shadow=True
)
plt.savefig("linbladian_spectrum.svg")#,dpi=500)
plt.show()




import numpy as np;
import scipy.signal as signal;
import time;    # For perf_counter to compare execution times of SLL-suppression algorithms

# To run only one specific (SNR, delta) scenario:
# 1. Remove the lines "for Pn, SNRValue in zip(Noise, SNR):" and "for d in delta:" and get rid of the respective indentations
# 2. Set the "SNR" and "delta" variables as your desired scalars
# 3. Rename "Noise" to "Pn", "delta" to "d" and "SNR" to "SNRValue"

# Initializations
M = 24;
N = 5;
bd = np.pi;
IdM = np.eye(M);

SNR = np.array([-10, 0, 10, 20]);
Ps = 1;
Noise = Ps/(10**(SNR/10));
Rgigi = Ps*np.eye(N);

delta = [6, 8, 10, 12, 14, 16];
AngleStart = 30;
AngleEnd = 150;
jbd = 1j*bd*np.arange(M)[:, np.newaxis];
theta = np.arange(0, 180, 0.1);     # 0.1 deg resolution
ath = np.exp(jbd @ np.cos(np.radians(theta[np.newaxis, :])));

# Creation of all required permutations of (0, 1, 2, 3, 4, 5) to be used with "d" to generate the angles
FirstAngleSet = np.arange(N+1);
Mask = FirstAngleSet[np.newaxis, :] != FirstAngleSet[:, np.newaxis];
Others = np.tile(FirstAngleSet, (N+1, 1))[Mask].reshape(N+1, N);
AnglePermutations = np.column_stack((FirstAngleSet, Others))[np.newaxis, :, :];

with open("results.csv", "w", encoding="utf-8") as ResultsFile:
    print("# MVDR", file=ResultsFile);
with open("iterations.csv", "w", encoding="utf-8") as IterationsFile:
    print("# MVDR", file=IterationsFile);

AoAdev_SINR_SLL = np.empty((1, 2*(N+1)+2));
DesiredSLL = -20;
Margin = 10**(-12);         # Used to exclude the main lobe from being considered a side-lobe when finding peaks in the array factor
LoadingFactor = 0.125;      # For diagonal loading of Rxx to prevent ill-conditioning errors. Value may not be optimal.
BatchSize = 2;              # Number of side-lobes with height above the desired SLL to attempt to squash in each iteration
Counter = 0;                # Counts the total number of times (in all scenarios, in all permutations, etc.) the SLL was above the desired SLL
TotalCounter = 0;           # Counts how many different cases are ran in total

Start = time.perf_counter();
for Pn, SNRValue in zip(Noise, SNR):
    Rnn = Pn*IdM;
    with open("results.csv", "a", encoding="utf-8") as ResultsFile:
        print(f'# SNR={SNRValue},Δθ0(min),Δθ0(max),Δθ0(mean),Δθ0(std),Δθ(min),Δθ(max),Δθ(mean),Δθ(std),SINR(min),SINR(max),SINR(mean),SINR(std),SLL(min),SLL(max),SLL(mean),SLL(std)',
              file=ResultsFile);
    with open("iterations.csv", "a", encoding="utf-8") as IterationsFile:
        print(f'# SNR={SNRValue},Iterations(min),Iterations(max),Iterations(mean),Iterations(std)', file=IterationsFile);
    for d in delta:
        # (1) Create all required angle combinations
        Angles = np.arange(AngleStart, AngleEnd-N*d+1)[:, np.newaxis, np.newaxis] + d*AnglePermutations;
        Iterations = np.zeros((Angles.shape[0], N+1), dtype=int);
        with open("AoAdev_SINR_SLL.txt", "w", encoding="utf-8") as AoAdev_SINR_SLL_File:
            print("# θ0,θ1,θ2,θ3,θ4,θ5,Δθ0,Δθ1,Δθ2,Δθ3,Δθ4,Δθ5,SINR,SLL", file=AoAdev_SINR_SLL_File);
        for th in range(Angles.shape[0]):
            for Permutation in range(N+1):
                TotalCounter += 1;
                A = np.exp(jbd @ np.cos(np.radians(Angles[th, Permutation, np.newaxis, :])));
                ad = A[:, 0, np.newaxis];
                Ai = A[:, 1:6];
                Ruu = Ai @ Rgigi @ Ai.conj().T + Rnn;
                Rdd = Ps * ad @ ad.conj().T;

                Rxx = A @ (Ps*np.eye(A.shape[1])) @ A.conj().T + Rnn;
                Rxx = Rxx + LoadingFactor*(np.trace(Rxx)/M)*IdM;    # Diagonal loading
                w = np.linalg.solve(Rxx, ad);
                AF = np.abs(w.conj().T @ ath).squeeze();
                AFdB = 20*np.log10(AF/AF.max());
                Peaks = signal.find_peaks(AFdB, height=(DesiredSLL + Margin, -Margin))[0];

                # (2) Fake null placement for SLL suppression
                # While there are still peaks above the desired SLL and the array can accommodate more interferers:
                # 1. Find them
                # 2. If there's at least 2 and the array can accommodate 2 more interferers, add the 2 largest as fake nulls
                # 3. If there's only 1 or the array can't accommodate 2 more interferers, add the largest one as a fake null
                # The idea of adding two at a time is to speed up the process, and also to perhaps take advantage of the fact
                # that the two largest side lobes often appear as a pair, one on each side of the main lobe, and thus squashing both 
                # at once may not be a terrible idea. Larger batch sizes do not appear to work so well.
                while len(Peaks) > 0 and A.shape[1] < M:
                    Iterations[th, Permutation] += 1;
                    # To test the one-at-a-time approach:
                    # 1. Comment the "if" and "else" lines and everything in between
                    # 2. Remove the indentation from the two lines that follow the "else" line
                    if len(Peaks) >= BatchSize and A.shape[1] + BatchSize <= M:
                        NewInterferers = theta[Peaks[AFdB[Peaks].argpartition(-BatchSize)[-BatchSize:]]];
                        A = np.column_stack((A, np.exp(jbd @ np.cos(np.radians(NewInterferers[np.newaxis, :])))));
                    else:
                        NewInterferer = theta[Peaks[AFdB[Peaks].argmax()]];
                        A = np.column_stack((A, np.exp(jbd * np.cos(np.radians(NewInterferer)))));
                    # (3) i. Calculate weights using MVDR algorithm
                    Rxx = A @ (Ps*np.eye(A.shape[1])) @ A.conj().T + Rnn;
                    Rxx = Rxx + LoadingFactor*(np.trace(Rxx)/M)*IdM;    # Diagonal loading
                    w = np.linalg.solve(Rxx, ad);
                    # (3) ii. Calculate array factor
                    AF = np.abs(w.conj().T @ ath).squeeze();
                    AFdB = 20*np.log10(AF/AF.max());
                    Peaks = signal.find_peaks(AFdB, height=(DesiredSLL + Margin, -Margin))[0];
                
                # (3) iii. b. Calculate SLL
                SLL = AFdB[signal.find_peaks(AFdB, height=(-100, -Margin))[0]].max();
                if SLL > DesiredSLL:
                    Counter += 1;
                # (3) ii. Calculate SINR. The calculation does not take into account fake interferers;
                # it is assumed that "SINR" only cares about real interferers and noise. Thus, the addition
                # of fake nulls does not change the Rdd or Ruu matrices, therefore the weights become the only differentiating factor.
                SINR = 10*np.log10(np.squeeze((w.conj().T @ Rdd @ w) / (w.conj().T @ Ruu @ w)));
                Zeroes = theta[signal.find_peaks(-AFdB)[0]];    # Find all zero positions for the upcoming AoA deviation calculations
                AoAdev_SINR_SLL[0, 0:6] = Angles[th, Permutation, 0:6];
                # (3) iii. a. Calculate AoA deviations
                AoAdev_SINR_SLL[0, 6] = np.abs(theta[AFdB.argmax()] - Angles[th, Permutation, 0]);                      # Main lobe deviation
                AoAdev_SINR_SLL[0, 7:12] = (np.abs(Zeroes[:, np.newaxis] - Angles[th, Permutation, 1:6])).min(axis=0);  # Null deviations
                AoAdev_SINR_SLL[0, 12:14] = [np.real(SINR), SLL];
                # (4) Save above to file named "AoAdev_SINR_SLL.txt"
                with open("AoAdev_SINR_SLL.txt", "a", encoding="utf-8") as AoAdev_SINR_SLL_File:
                    np.savetxt(AoAdev_SINR_SLL_File, AoAdev_SINR_SLL, delimiter=",", fmt='%g', encoding="utf-8");
        
        # (5) Load file and calculate required metrics
        Results = np.loadtxt("AoAdev_SINR_SLL.txt", delimiter=",", encoding="utf-8");
        Results = np.array([d, Results[:, 6].min(), Results[:, 6].max(), Results[:, 6].mean(), Results[:, 6].std(),
                            Results[:, 7:12].min(), Results[:, 7:12].max(), Results[:, 7:12].mean(), Results[:, 7:12].std(),
                            Results[:, 12].min(), Results[:, 12].max(), Results[:, 12].mean(), Results[:, 12].std(),
                            Results[:, 13].min(), Results[:, 13].max(), Results[:, 13].mean(), Results[:, 13].std()])[np.newaxis, :];
        # Saving to file helps with experimentation and transferring to report
        with open("results.csv", "a", encoding="utf-8") as ResultsFile:
            print("δ=", end='', file=ResultsFile);
            np.savetxt(ResultsFile, Results, delimiter=",", fmt='%.3f', encoding="utf-8");
        # Extras
        Iterations = np.array([d, Iterations.min(), Iterations.max(), Iterations.mean(), Iterations.std()])[np.newaxis, :];
        with open("iterations.csv", "a", encoding="utf-8") as IterationsFile:
            print("δ=", end='', file=IterationsFile);
            np.savetxt(IterationsFile, Iterations, delimiter=",", fmt='%.3f', encoding="utf-8");
End = time.perf_counter();

print(f'The SLL was above {DesiredSLL} dB in {Counter} out of {TotalCounter} scenarios, or {(Counter/TotalCounter)*100:.2f}% of the time.');
print(f'Time taken for program execution: {End-Start:.2f} seconds');

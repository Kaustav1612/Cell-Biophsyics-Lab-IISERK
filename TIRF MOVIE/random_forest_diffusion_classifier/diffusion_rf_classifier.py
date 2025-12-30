import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from sklearn.base import BaseEstimator

# ==============================================
# Part 1: Synthetic Data Generation
# ==============================================
def generate_trajectory(alpha, n_frames=100, D=1.0):
    displacements = np.random.normal(0, np.sqrt(2*D), (n_frames-1, 2))
    if alpha != 1:
        time_scaling = np.arange(1, n_frames)**((alpha-1)/2)
        displacements = displacements * time_scaling[:, np.newaxis]
    return np.vstack(([0, 0], np.cumsum(displacements, axis=0)))

def calculate_ta_msd(trajectory, max_lag=None):
    n_frames = len(trajectory)
    max_lag = max_lag or n_frames//4
    msd = np.zeros(max_lag)
    for lag in range(1, max_lag+1):
        diffs = trajectory[lag:] - trajectory[:-lag]
        msd[lag-1] = np.mean(np.sum(diffs**2, axis=1))
    return np.arange(1, max_lag+1), msd

def extract_features(trajectory):
    taus, msd = calculate_ta_msd(trajectory)
    N = len(trajectory)
    features = []
    
    # MSD fitting features
    try:
        popt, _ = curve_fit(lambda t, a, b: b*t**a, taus, msd, maxfev=1000)
        features.extend([popt[0], np.log(popt[1]) if popt[1] > 0 else -10])
    except:
        features.extend([np.nan, np.nan])
    
    # Maximum distance features
    D_N = np.max(np.linalg.norm(trajectory - trajectory[0], axis=1))
    step_sizes = np.linalg.norm(np.diff(trajectory, axis=0), axis=1)
    σ_N_sq = np.sum(step_sizes**2) / (2*N)
    features.append(D_N / np.sqrt(σ_N_sq * (N-1)) if σ_N_sq > 0 else np.nan)
    
    # p-variation features
    for p in range(1, 6):
        n = max(1, N//10)
        V_p = sum(np.linalg.norm(trajectory[(k+1)*n] - trajectory[k*n])**p 
                 for k in range(N//n-1))
        features.append(V_p)
    
    # Additional MSD features
    features.extend([
        np.log(msd[0]) if len(msd) > 0 and msd[0] > 0 else -10,
        np.var(msd) if len(msd) > 1 else 0
    ])
    
    return features

# ==============================================
# Part 2: Data Preparation and Training
# ==============================================
alphas = [0.3, 0.7, 1.0, 1.3]
class_labels = ['subdiffusive', 'confined', 'normal', 'superdiffusive']  # Discrete labels
traj_num =1000
X, y = [], []
for _ in range(traj_num):
    index  = np.random.randint(0,3)
    traj = generate_trajectory(alphas[index])
    features = extract_features(traj)
    if not np.isnan(features).any():  
        X.append(features)
        y.append(class_labels[index])

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize and train classifier
clf = RandomForestClassifier(
    n_estimators=200,
    class_weight='balanced',
    random_state=42,
    n_jobs=-1
)
clf.fit(X_train, y_train)

# Evaluation
print("Classification Performance:")
print(classification_report(y_test, clf.predict(X_test)))




# Load the Excel data
df = pd.read_excel('trajectories_data_random.xlsx')

# Function to extract trajectories from the DataFrame
def extract_trajectories(df):
    trajectories = []
    particle_ids = set()
    
    # Get unique particle IDs from column names
    for col in df.columns:
        if '_X' in col:
            particle_id = col.split('_')[0][1:]  # Extract P1 from P1_X
            particle_ids.add(particle_id)
    
    # For each particle, combine X, Y, and Frame data
    for pid in sorted(particle_ids):
        x_col = f'P{pid}_X'
        y_col = f'P{pid}_Y'
        frame_col = f'P{pid}_Frame'
        
        if x_col in df.columns and y_col in df.columns:
            # Combine valid positions (drop NaN rows)
            particle_data = df[[frame_col, x_col, y_col]].dropna()
            if len(particle_data) > 1:  # Need at least 2 points for MSD
                # Convert to Nx2 array (X,Y positions)
                positions = particle_data[[x_col, y_col]].values
                trajectories.append(positions)
    
    return trajectories

# MSD calculation function (needed for your classifier)
def calculate_ta_msd(traj):
    """Calculate time lags and MSD for a trajectory"""
    n = len(traj)
    taus = np.arange(1, min(n, 50))  
    msd = np.zeros_like(taus, dtype=float)
    
    for i, tau in enumerate(taus):
        if tau < n:
            disp = traj[tau:] - traj[:-tau]
            msd[i] = np.mean(disp[:,0]**2 + disp[:,1]**2)
    
    return taus, msd

def classify_real_trajectories(trajectories, clf):
    """Classify experimental trajectories using the full feature set"""
    all_features = []
    valid_indices = [] # Keep track of which trajectories yielded valid features

    for i, traj in enumerate(trajectories):
        current_features = extract_features(traj) # Call your existing feature extraction
        
        # Check for NaN values before appending
        if not np.isnan(current_features).any():  
            all_features.append(current_features)
            valid_indices.append(i) # Store index of valid trajectory
        # else:
        #   You could also append NaNs if your classifier is set up to handle them,
        #   or simply skip the invalid trajectory as done before.
        #   For now, we'll continue to filter them out.
            
    if not all_features: # Handle case where no valid trajectories are found
        print("Warning: No valid features extracted from real trajectories.")
        return np.array([]), np.array([])
        
    # Convert list of lists to a NumPy array for prediction
    feature_array = np.array(all_features)
    
    # Predict the classes
    predictions = clf.predict(feature_array)
    
    # Return predictions and the indices of the trajectories that were actually classified
    return predictions, np.array(valid_indices)

# Visualization function
def plot_msd_classification(trajectories, alphas):
    plt.figure(figsize=(12, 6))
    colors = plt.cm.viridis(np.linspace(0, 1, len(np.unique(alphas))))
    
    for alpha, c in zip(np.unique(alphas), colors):
        idx = np.where(alphas == alpha)[0]
        for i in idx[:10]:  # Plot first 10 examples per class
            taus, msd = calculate_ta_msd(trajectories[i])
            plt.loglog(taus, msd, color=c, alpha=0.3, 
                      label=f'α={alpha}' if i == idx[0] else None)
    
    plt.xlabel('Time lag (τ)')
    plt.ylabel('MSD (nm²)')
    plt.title('TA-MSD Colored by Predicted α')
    plt.legend()
    plt.show()

# Main workflow
if __name__ == "__main__":
    # ... (previous code remains the same) ...

    # 1. Load and prepare trajectories
    real_trajectories = extract_trajectories(df)
    print(f"Number of real trajectories extracted: {len(real_trajectories)}")

    # 3. Classify trajectories
    # The 'valid' array will now contain indices corresponding to original real_trajectories
    predicted_alphas, valid_original_indices = classify_real_trajectories(real_trajectories, clf)
    
    print(f"Number of trajectories classified: {len(predicted_alphas)}")

    # 4. Visualize results
    # Pass only the trajectories that were successfully classified
    classified_trajectories_for_plot = [real_trajectories[i] for i in valid_original_indices]
    
    # Ensure predicted_alphas is passed correctly for plotting (it's the classes)
    plot_msd_classification(classified_trajectories_for_plot, predicted_alphas)





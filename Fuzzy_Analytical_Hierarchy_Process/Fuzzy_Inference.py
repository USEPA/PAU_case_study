# def _fuzzy_inference(self, df):
#     N = dict()
#     rows = df.shape[0]
#     if rows == 1:
#         return df
#     else:
#         naics = df['PRIMARY NAICS CODE'].unique().tolist()[0]
#         cas = df['CAS NUMBER'].unique().tolist()[0]
#         efficiency = df['EFFICIENCY RANGE CODE'].tolist()
#         concentration = df['RANGE INFLUENT CONCENTRATION'].tolist()
#         waste_code = df['WASTE STREAM CODE'].tolist()
#         n_actual = len(waste_code)
#         if self.Year <= 2004:
#             efficiency = df['EFFICIENCY ESTIMATION'].tolist()
#             #indices = [i if x == 0.0 else None for i, x in enumerate(df['EFFICIENCY ESTIMATION'].tolist())]
#             indices = [i if x == 0.0 else None for i, x in enumerate(efficiency)]
#             indices_f = {v: i for i, v in enumerate(set(range(len(waste_code))) - set(indices))}
#             efficiency = [x for i, x in enumerate(efficiency) if not i in indices]
#             concentration = [x for i, x in enumerate(concentration) if not i in indices]
#             waste_code = [x for i, x in enumerate(waste_code) if not i in indices]
#         n = len(waste_code)
#         if n >= 2:
#             # Based on position
#             m_criteria = 1
#             N_aux = np.empty((n, n))
#             N_aux[:] = np.nan
#             for i in range(n):
#                 for j in range(i, n):
#                     dist = j - i
#                     if dist == 0:
#                         N_aux[i][j] = 0
#                     elif dist == 4:
#                         N_aux[i][j] = 4
#                     elif dist == 3:
#                         N_aux[i][j] = 3
#                     elif dist == 2:
#                         N_aux[i][j] = 2
#                     elif dist == 1:
#                         N_aux[i][j] = 1
#             N =  N_aux
#             # Based on efficiency
#             m_criteria = m_criteria + 1
#             #dict_eff = {'E1':0, 'E2':0, 'E3':1, 'E4':2, 'E5':3, 'E6':4}
#             N_aux = np.empty((n, n))
#             N_aux[:] = np.nan
#             for i in range(n):
#                 for j in range(i, n):
#                     #dist = dict_eff[efficiency[j]] - dict_eff[efficiency[i]]
#                     if self.Year >= 2005:
#                         dist =  (self._range_to_efficiency_estimation(efficiency[i]) - \
#                                 self._range_to_efficiency_estimation(efficiency[j]))/100
#                     else:
#                         dist = (efficiency[i] - efficiency[j])/100
#                     if dist == 0:
#                         N_aux[i][j] = 0
#                     elif (0.5 < dist) & (dist <= 1):#dist == 4:
#                         N_aux[i][j] = 4
#                     elif (0.245 < dist) & (dist <= 0.5):#dist == 3:
#                         N_aux[i][j] = 3
#                     elif (0.003 < dist) & (dist <= 0.245):#dist == 2:
#                         N_aux[i][j] = 2
#                     elif (0 < dist) & (dist <= 0.003):#dist == 1:
#                         N_aux[i][j] = 1
#                     elif (-1 <= dist) & (dist < -0.5):#dist == -4:
#                         N_aux[i][j] = -4
#                     elif (-0.5 <= dist) & (dist < -0.245):#dist == -3:
#                         N_aux[i][j] = -3
#                     elif (-0.245 <= dist) & (dist < -0.003):#dist == -2:
#                         N_aux[i][j] = -2
#                     elif (-0.003 <= dist) & (dist < 0):#dist == -1:
#                         N_aux[i][j] = -1
#             N =  np.concatenate((N, N_aux), axis = 0)
#             # Based on concentration
#             if self.Year <= 2004:
#                 m_criteria = m_criteria + 1
#                 #dict_con = {1:0, 2:1, 3:2, 4:3, 5:4}
#                 N_aux = np.empty((n, n))
#                 N_aux[:] = np.nan
#                 for i in range(n):
#                     for j in range(i, n):
#                         dist =  (self._range_to_concentration_estimation(concentration[i]) - \
#                                 self._range_to_concentration_estimation(concentration[j]))/100
#                         #dist = dict_con[concentration[j]] - dict_con[concentration[i]]
#                         if dist == 0:
#                             N_aux[i][j] = 0
#                         elif (0.25 < dist) & (dist <= 0.5):#dist == 4:
#                             N_aux[i][j] = 4
#                         elif (0.125 < dist) & (dist <= 0.25):#dist == 3:
#                             N_aux[i][j] = 3
#                         elif (0.004 < dist) & (dist <= 0.125):#dist == 2:
#                             N_aux[i][j] = 2
#                         elif (0 < dist) & (dist <= 0.004):#dist == 1:
#                             N_aux[i][j] = 1
#                         elif (-0.5 <= dist) & (dist < -0.25):#dist == -4:
#                             N_aux[i][j] = -4
#                         elif (-0.25 <= dist) & (dist < -0.125):#dist == -3:
#                             N_aux[i][j] = -3
#                         elif (-0.125 <= dist) & (dist < -0.004):#dist == -2:
#                             N_aux[i][j] = -2
#                         elif (-0.004 <= dist) & (dist < 0):#dist == -1:
#                             N_aux[i][j] = -1
#                 N =  np.concatenate((N, N_aux), axis = 0)
#             # Based on statistics
#             m_criteria = m_criteria + 1
#             N_aux =  np.empty((n, n))
#             N_aux[:] = np.nan
#             df_statistics = pd.read_csv(self._dir_path + '/Statistics/DB_for_Statistics_' + str(self.Year) + '.csv',
#                             usecols = ['CAS', 'NAICS', 'WASTE'],
#                             low_memory = False,
#                             converters = {'NAICS':  lambda x: str(int(float(x)))})
#             df_statistics = df_statistics.loc[df_statistics['CAS'] == cas]
#             if df_statistics.shape[0] != n:
#                 df_statistics['NAICS STRUCTURE'] = df_statistics.apply(lambda x: self._searching_naics(\
#                                     x['NAICS'], naics), axis = 1)
#                 naics_structure = ['National Industry', 'NAICS Industry', 'Industry Group',
#                                  'Subsector', 'Sector', 'Nothing']
#                 count = {}
#                 j = 0
#                 while not count:
#                     structure = naics_structure[j]
#                     df_naics = df_statistics.loc[df_statistics['NAICS STRUCTURE'] ==  structure]
#                     if (structure == 'National Industry') and (df_naics.shape[0] == n):
#                         count = {i: waste_code.count(i)/len(waste_code) for i in list(set(waste_code))}
#                     else:
#                         values = df_naics['WASTE'].value_counts(normalize = True).keys().tolist()
#                         counts = df_naics['WASTE'].value_counts(normalize = True).tolist()
#                         count = {values[i]: counts[i] for i in range(len(values))}
#                     j = j + 1
#             else:
#                 count = {i: waste_code.count(i)/len(waste_code) for i in list(set(waste_code))}
#             for i in range(n):
#                 for j in range(i, n):
#                     dist =  count[waste_code[i]] - count[waste_code[j]]
#                     if dist == 0.0:
#                         N_aux[i][j] = 0
#                     elif (dist > 0.0) and (dist <= 0.25):
#                         N_aux[i][j] = 1
#                     elif (dist > 0.25) and (dist <= 0.5):
#                         N_aux[i][j] = 2
#                     elif (dist > 0.5) and (dist <= 0.75):
#                         N_aux[i][j] = 3
#                     elif (dist > 0.75) and (dist <= 1.0):
#                         N_aux[i][j] = 4
#                     elif (dist < 0.0) and (dist >= -0.25):
#                         N_aux[i][j] = -1
#                     elif (dist < -0.25) and (dist >= -0.5):
#                         N_aux[i][j] = -2
#                     elif (dist < -0.5) and (dist >= -0.75):
#                         N_aux[i][j] = -3
#                     elif (dist < -0.75) and (dist >= -1.0):
#                         N_aux[i][j] = -4
#             N =  np.concatenate((N, N_aux), axis = 0)
#             # Definition of variables
#             W = np.zeros((1, n)) # initial value of weights vector (desired output)
#             w = np.zeros((m_criteria, n)) # initial value of weithts matrix
#             phi = np.zeros((m_criteria, 1)) # diversification degree
#             eps = np.zeros((m_criteria, 1)) # entropy
#             theta = np.zeros((1, m_criteria)) # criteria' uncertainty degrees
#             sumphi = 0 # Initial value of diversification degree
#             # Triangular fuzzy numbers (TFN)
#             # TFN is a vector with n segments and each one has 3 different numbers
#             # The segments are the linguistic scales
#             # The 3 differente numbers in the segments are a triangular fuzzy number (l,m,u)
#             # the first segment: equal importance; the second one: moderate importance of one over another;
#             # the third one: strong importance of one over another; the fourth one: very strong importance of one over another
#             # the fifth one: Absolute importance of one over another
#             TFN = np.array([1, 1, 1 ,2/3, 1, 3/2, 3/2, 2, 5/2, 5/2, 3, 7/2, 7/2, 4, 9/2])
#             for k in range(1, m_criteria + 1):
#                 a = np.zeros((1, n*n*3))  # Comparison matrix (In this case is a vector because of computational memory
#                 for i in range(k*n-(n-1), k*n + 1):
#                     for j in range(i-n*(k-1), n + 1):
#                         # This is the position of the third element of the segment for
# 			            # a*(i,j) (upper triangular part)
#                         jj = 3*(n*((i-n*(k-1))-1)+j)
#                         # This is the position of the thrid element of the segment for
# 			            # a*(j,i) (lower triangular part)
#                         jjj = 3*(n*(j-1) + i-n*(k-1))
#                         if N[i - 1][j - 1] == -4:
#                             a[0][jjj-3:jjj] =  TFN[12:15]
#                             a[0][jj-3:jj] = np.array([TFN[14]**-1, TFN[13]**-1, TFN[12]**-1])
#                         elif N[i - 1][j - 1] == -3:
#                             a[0][jjj-3:jjj] =  TFN[9:12]
#                             a[0][jj-3:jj] = np.array([TFN[11]**-1, TFN[10]**-1, TFN[9]**-1])
#                         elif N[i - 1][j - 1] == -2:
#                             a[0][jjj-3:jjj] = TFN[6:9]
#                             a[0][jj-3:jj] = np.array([TFN[8]**-1, TFN[7]**-1, TFN[6]**-1])
#                         elif N[i - 1][j - 1] == -1:
#                             a[0][jjj-3:jjj] =  TFN[3:6]
#                             a[0][jj-3:jj] = np.array([TFN[5]**-1, TFN[4]**-1, TFN[3]**-1])
#                         elif N[i - 1][j - 1] == 0:
#                             a[0][jj-3:jj] =  TFN[0:3]
#                             a[0][jjj-3:jjj] = TFN[0:3]
#                         elif N[i - 1][j - 1] == 1:
#                             a[0][jj-3:jj] =  TFN[3:6]
#                             a[0][jjj-3:jjj] = np.array([TFN[5]**-1, TFN[4]**-1, TFN[3]**-1])
#                         elif N[i - 1][j - 1] == 2:
#                             a[0][jj-3:jj] = TFN[6:9]
#                             a[0][jjj-3:jjj] = np.array([TFN[8]**-1, TFN[7]**-1, TFN[6]**-1])
#                         elif N[i - 1][j - 1] == 3:
#                             a[0][jj-3:jj] =  TFN[9:12]
#                             a[0][jjj-3:jjj] = np.array([TFN[11]**-1, TFN[10]**-1, TFN[9]**-1])
#                         elif N[i - 1][j - 1] == 4:
#                             a[0][jj-3:jj] =  TFN[12:15]
#                             a[0][jjj-3:jjj] = np.array([TFN[14]**-1, TFN[13]**-1, TFN[12]**-1])
#                 # (2) fuzzy synthetic extension
#                 A = np.zeros((n,3))
#                 B = np.zeros((1,3))
#                 for i in range(1, n + 1):
#                     for j in range(1, n + 1):
#                         jj = 3*(n*(i-1)+j)
#                         A[i - 1][:] = A[i - 1][:] + a[0][jj-3:jj]
#                     B = B + A[i - 1][:]
#                 BB = np.array([B[0][2]**-1, B[0][1]**-1, B[0][0]**-1])
#                 S = A*BB
#                 # (3) Degree of possibility
#                 for i in range(n):
#                     V = np.zeros(n)
#                     for j in range(n):
#                         if S[i][1] >= S[j][1]:
#                             V[j] = 1
#                         elif S[j][0] >= S[i][2]:
#                             V[j] = 0
#                         else:
#                             V[j] = (S[j][0] - S[i][2])/((S[i][1] - S[i][2]) - (S[j][1] - S[j][0]))
#                     w[k - 1][i] = np.min(V)
#                 # (4) Weight of each flow for a criterium
#                 w[k - 1][:] = (np.sum(w[k - 1][:])**-1)*w[k - 1][:]
#                 # (5) Criteria' uncertainty degrees
#                 for i in range(n):
#                     if w[k - 1][i] != 0:
#                         eps[k - 1][0] = eps[k - 1][0] - w[k - 1][i]*np.log(w[k - 1][i])
#                 eps[k - 1][0] = eps[k - 1][0]/np.log(n)
#                 phi[k - 1][0] = 1 + eps[k - 1]
#                 sumphi = sumphi + phi[k - 1][0]
#             # (6) Final weight of all flows
#             for i in range(n):
#                 for k in range(m_criteria):
#                     theta[0][k] = phi[k][0]/sumphi
#                     W[0][i] = W[0][i] + w[k][i]*theta[0][k]
#             W = (np.sum(W)**-1)*W # Final weights
#             if self.Year <= 2004:
#                 W = np.array([0.0 if i in indices else W[0][indices_f[i]] for i in range(n_actual)])
#         else:
#             W = np.array([0.0 if i in indices else 1.0 for i in range(n_actual)])
#         df = df.assign(WEIGHT = W)
#         df['MANAGED FLOW'] = df.apply(lambda x: x['MANAGED FLOW']*x['WEIGHT'], axis = 1)
#         df.drop(columns = 'WEIGHT', inplace = True)
#         return df

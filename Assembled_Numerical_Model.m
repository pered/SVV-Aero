load('Cit_par.mat')
import Symmetrical_Model_Numerical

xu = V0 * CXu / (c * 2 * muc);
xa = V0 * CXa / (c * 2 * muc);
x0 = V0 * CZ0 / (c * 2 * muc);

zu = V0 * CZu / (c * (2*muc - CZadot));
za = V0 * CZa / (c * (2*muc - CZadot));
z0 = V0 * CZ0 / (c * (2*muc - CZadot));
zq = V0 * CZq / (c * (2*muc - CZadot));

mu = V0 * (Cmu + CZu * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2);
ma = V0 * (Cma + CZa * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2);
m0 = V0 * (CX0 * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2);
mq = V0 * (Cmq + CZq * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2);

Citation = Symmetrical_Model_Numerical;
/* Copyright (C) 2013 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Copyright (c) 2010, Moeti Ncube, All rights reserved. (original code)
 * (http://www.mathworks.com/matlabcentral/fileexchange/29751-calibration-method-for-the-schwartz-smith-model)
 * Author(s): Stefano Angeleri (stefano.angeleri@fairmat.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Original code from moeti ncube is licensed under the BSD license:
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the distribution
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

using System;
using DVPLI;
using Fairmat.Math.Algebra;
using Range = DVPLI.Range;

namespace SchwartzSmithModel
{
    class MainClass
    {
        // TODO: fe[i,2] =   NORM is missing!!!
        public static Matrix CM(Matrix mu, Matrix sumP21, Matrix sumP1)
        {
            /*
             * function y=cm(mu,sumP21,sumP1)
                mu(end+1,:)=ones(1,length(mu));
                mu1=mu(:,1:end-1);
                mu2=mu(:,2:end);
                sumP21(end+1,1:end)=0;
                sumP21(1:end,end+1)=0;
                sumP1(end+1,1:end)=0;
                sumP1(1:end,end+1)=0;
                Ahat=(sumP21+mu2*mu1')*inv((sumP1+mu1*mu1'));
                y=Ahat(1:end-1,end);
            */
            Matrix muTmp = mu.Clone();
            muTmp[muTmp.R, Range.All] = Matrix.Ones(1, muTmp.C);
            Matrix mu1 = muTmp[Range.All, new Range(0, muTmp.C - 2)];
            Matrix mu2 = muTmp[Range.All, new Range(1, muTmp.C - 1)];
            Matrix sumP21Tmp = sumP21.Clone();
            Matrix sumP1Tmp = sumP1.Clone();
            sumP21Tmp.Resize(sumP21Tmp.R + 1, sumP21Tmp.C + 1);
            sumP1Tmp.Resize(sumP1Tmp.R + 1, sumP1Tmp.C + 1);
            Matrix Ahat = (sumP21Tmp + mu2 * mu1.Transpose()) * (sumP1Tmp + mu1 * mu1.Transpose()).Inverse();
            return Ahat[new Range(0, Ahat.R - 2), Ahat.C - 1];
        }

        public static Matrix HM(Matrix x, Matrix mu, Matrix sumP)
        {
            /*
             * function y=hm(x,mu,sumP)
                x(end+1,:)=ones(1,length(mu));
                mu(end+1,:)=ones(1,length(mu));
                sumP(end+1,1:end)=0;
                sumP(1:end,end+1)=0;
                Hhat=(x*mu')*inv(sumP+mu*mu');
                y=Hhat(1:end-1,1:end-1);
            */
            Matrix xTmp = x.Clone();
            Matrix muTmp = mu.Clone();
            Matrix sumPTmp = sumP.Clone();
            xTmp[xTmp.R, Range.All] = Matrix.Ones(1, muTmp.C);
            muTmp[muTmp.R, Range.All] = Matrix.Ones(1, muTmp.C);
            sumPTmp.Resize(sumPTmp.R + 1, sumPTmp.C + 1);
            Matrix Hhat = (xTmp * muTmp.Transpose()) * ((sumPTmp + muTmp * muTmp.Transpose()).Inverse());
            return Hhat[new Range(0, Hhat.R - 2), new Range(0, Hhat.C - 2)];
        }

        public static Matrix AM(Matrix mu, Matrix sumP21, Matrix sumP1)
        {
            /*function y=am(mu,sumP21,sumP1)
                mu(end+1,:)=ones(1,length(mu));
                mu1=mu(:,1:end-1);
                mu2=mu(:,2:end);
                sumP21(end+1,1:end)=0;
                sumP21(1:end,end+1)=0;
                sumP1(end+1,1:end)=0;
                sumP1(1:end,end+1)=0;
                Ahat=(sumP21+mu2*mu1')*inv((sumP1+mu1*mu1'));
                y=Ahat(1:end-1,1:end-1);
            */
            Matrix muTmp = mu.Clone();
            muTmp[muTmp.R, Range.All] = Matrix.Ones(1, muTmp.C);
            Matrix mu1 = muTmp[Range.All, new Range(0, muTmp.C - 2)];
            Matrix mu2 = muTmp[Range.All, new Range(1, muTmp.C - 1)];
            Matrix sumP21Tmp = sumP21.Clone();
            Matrix sumP1Tmp = sumP1.Clone();
            sumP21Tmp.Resize(sumP21Tmp.R + 1, sumP21Tmp.C + 1);
            sumP1Tmp.Resize(sumP1Tmp.R + 1, sumP1Tmp.C + 1);
            Matrix Ahat = (sumP21Tmp + mu2 * mu1.Transpose()) * (sumP1Tmp + mu1 * mu1.Transpose()).Inverse();
            return Ahat[new Range(0, Ahat.R - 2), new Range(0, Ahat.C - 2)];
        }

        public static Matrix QM(Matrix x, Matrix dhat, Matrix Hhat, Matrix mu, Matrix sumP, int T)
        {
            /*x(end+1,:)=ones(1,length(mu));
                mu(end+1,:)=ones(1,length(mu));
                sumP(end+1,1:end)=0;
                sumP(1:end,end+1)=0;
                %Hhat=[Hhat,dhat;0,0,1];
                Hhat=[Hhat,dhat;zeros(1,size(Hhat,2)),1];
                Qhat=((x-Hhat*mu)*(x-Hhat*mu)'+Hhat*sumP*Hhat')/T;
                y=Qhat(1:end-1,1:end-1);
                */
            Matrix xTmp = x.Clone();
            Matrix sumPTmp = sumP.Clone();
            Matrix muTmp = mu.Clone();
            xTmp[xTmp.R, Range.All] = Matrix.Ones(1, muTmp.C);
            muTmp[muTmp.R, Range.All] = Matrix.Ones(1, muTmp.C);
            sumPTmp.Resize(sumPTmp.R + 1, sumPTmp.C + 1);
            Hhat = new Matrix(new Matrix[][] { new[] { Hhat, dhat }, new[] { new Matrix(1, Hhat.C, true), new Matrix(new double[,] { { 1 } }) } });
            Matrix Qhat = ((xTmp - Hhat * muTmp) * ((xTmp - Hhat * muTmp).Transpose()) + Hhat * sumPTmp * Hhat.Transpose()) / T;

            return Qhat[new Range(0, Qhat.R - 2), new Range(0, Qhat.C - 2)];
        }

        public static Matrix DM(Matrix obs, Matrix mu, Matrix sumP)
        {
            /*function y=dm(obs,mu,sumP)
            obs(end+1,:)=ones(1,length(mu));
            mu(end+1,:)=ones(1,length(mu));
            sumP(end+1,1:end)=0;
            sumP(1:end,end+1)=0;
            Hhat=(obs*mu')*inv(sumP+mu*mu');
            y=Hhat(1:end-1,end);*/
            Matrix obsTmp = obs.Clone();
            Matrix muTmp = mu.Clone();
            Matrix sumPTmp = sumP.Clone();
            obsTmp[obsTmp.R, Range.All] = Matrix.Ones(1, muTmp.C);
            muTmp[muTmp.R, Range.All] = Matrix.Ones(1, muTmp.C);
            sumPTmp.Resize(sumPTmp.R + 1, sumPTmp.C + 1);
            Matrix Hhat = (obsTmp * muTmp.Transpose()) * ((sumPTmp + muTmp * muTmp.Transpose()).Inverse());
            return Hhat[new Range(0, Hhat.R - 2), Hhat.C - 1];
        }

        public static Matrix WM(Matrix mu, Matrix sumP21, Matrix sumP2, Matrix sumP1, Matrix Ahat, Matrix Chat, int T)
        {
            /*
             * function y=wm(mu,sumP21,sumP2,sumP1,Ahat,Chat,T)
                mu(end+1,:)=ones(1,length(mu));
                mu1=mu(:,1:end-1);
                mu2=mu(:,2:end);
                sumP21(end+1,1:end)=0;
                sumP21(1:end,end+1)=0;
                sumP2(end+1,1:end)=0;
                sumP2(1:end,end+1)=0;
                sumP1(end+1,1:end)=0;
                sumP1(1:end,end+1)=0;
                %Ahat=[Ahat,Chat;0,0,1];
                Ahat=[Ahat,Chat;zeros(1,size(Ahat,2)),1];
                What=((mu2-Ahat*mu1)*(mu2-Ahat*mu1)'+sumP2-2*sumP21*Ahat'+Ahat*sumP1*Ahat')/(T-1);
                y=What(1:end-1,1:end-1);
                */
            Matrix muTmp = mu.Clone();
            muTmp[muTmp.R, Range.All] = Matrix.Ones(1, muTmp.C);
            Matrix sumP21Tmp = sumP21.Clone();
            Matrix sumP2Tmp = sumP2.Clone();
            Matrix sumP1Tmp = sumP1.Clone();
            Matrix mu1 = muTmp[Range.All, new Range(0, muTmp.C - 2)];
            Matrix mu2 = muTmp[Range.All, new Range(1, muTmp.C - 1)];
            sumP21Tmp.Resize(sumP21Tmp.R + 1, sumP21Tmp.C + 1);
            sumP2Tmp.Resize(sumP2Tmp.R + 1, sumP2Tmp.C + 1);
            sumP1Tmp.Resize(sumP1Tmp.R + 1, sumP1Tmp.C + 1);

            Ahat = new Matrix(new Matrix[][] { new[] { Ahat, Chat }, new[] { new Matrix(1, Ahat.C), new Matrix(new double[,] { { 1 } }) } });
            Matrix What = ((mu2 - Ahat * mu1) * ((mu2 - Ahat * mu1).Transpose()) + sumP2Tmp - 2 * sumP21Tmp * Ahat.Transpose() + Ahat * sumP1Tmp * Ahat.Transpose()) / (T - 1);

            // TODO: use a resize?
            return What[new Range(0, What.R - 2), new Range(0, What.C - 2)];
        }


        public static double TestOriginalSS(double dt, Matrix y, Vector maturities, Matrix x, Matrix x0, Matrix W0)
        {
            /*
             * function l = TestOriginalSS(dt,y,matur,x,x0,W0)

                % We rerun the estimate parameters given in the original schwartz-smith paper
                %(k,mue,sigmax,sigmae,pxe,lambdae,lambdax)

                theta0=[1.49,.03,.286,.145,.3,.0415,.157];
                % se0=[.03,.01,.144,.0728,.005,.0013,.044];
                s=[.042;.006;.003;.0000;.004];

                [C,A,W,d,H,Q]=matrices(theta0(1),theta0(2),theta0(3),theta0(4),theta0(5),...
                    theta0(6),theta0(7),s,matur,dt);

                %The Kalman Filter
                % x0=[.119;2.857]; %Initial state vector m(t)=E[xt;et]
                % W0=[.1,.1;.1,.1]; %Inital covariance matrix C(t)=cov[xt,et]
                % R=eye(size(W,1));
                [y_cond,v,a,a_cond,P,P_cond,F,logl]=kalman_filter(y',H,d,A,C,x0,W0,Q,W,0);
                l=sum(logl);
                est1=exp(a(:,1)+a(:,2));
                obs1=x(1:end-1);

                %Estimated MSE of SS-Opt estimates to observed data
                fnorm=norm(est1-obs1,'fro')/norm(obs1,'fro');
             */

            Vector theta0 = new Vector(new double[] { 1.49, .03, .286, .145, .3, .0415, .157 });
            Matrix s = new Matrix(new double[,] { { .042 }, { .006 }, { .003 }, { .0000 }, { .004 } });
            Matrix[] matOutput = Matrices(theta0[0], theta0[1], theta0[2], theta0[3], theta0[4], theta0[5], theta0[6], s, maturities, dt);
            Matrix C = matOutput[0];
            Matrix A = matOutput[1];
            Matrix W = matOutput[2];
            Matrix d = matOutput[3];
            Matrix H = matOutput[4];
            Matrix Q = matOutput[5];
            matOutput = KalmanFilter(y.Transpose(), H, d, A, C, x0, W0, Q, W, 0);
            Matrix y_cond = matOutput[0];
            Matrix v = matOutput[1];
            Matrix a = matOutput[2];
            Matrix a_cond = matOutput[3];
            Matrix P = matOutput[4];
            Matrix P_cond = matOutput[5];
            Matrix F = matOutput[6];
            Matrix logl = matOutput[7];
            double l = logl[Range.All,0].Sum();

            // TODO: This seems lost?
            Matrix est1 = Matrix.Exp(a[Range.All, 0] + a[Range.All, 1]);
            double obs1 = x[0, x.C - 2];
            // fnorm=norm(est1-obs1,'fro')/norm(obs1,'fro');
            return l;
        }

        public static Matrix Vecr(Matrix x)
        {
            Matrix xt = x.Transpose();
            Matrix output = new Matrix(x.R * x.C, 1);
            for (int c = 0, i = 0; c < xt.C; c++)
            {
                for (int r = 0; r < xt.R; r++)
                {
                    output[i, 0] = xt[r, c];
                    i++;
                }
            }

            return output;
        }

        public static Matrix[] KalmanSmoothing(Matrix y, Matrix Z, Matrix d, Matrix T, Matrix c, 
                                               Matrix a0, Matrix P0, Matrix H, Matrix Q, int timeVar)
        {
            /*
				function [a_smooth,P_smooth,P2_smooth] = Kalman_Smoothing(y,Z,d,T,c,a0,P0,H,Q,timevar)
				% Kalman smoothing
				% code is a Matlab translation of Thierry Roncalli's.

				nobs = size(y,1); 
				n    = size(y,2); 
				at   = a0; 
				Pt   = P0; 


				if timevar == 1
				    m=size(Z,2)/n; 
				    g=size(Q,2)/m;
				else
				    m=size(Z,2); 
				    g=size(Q,2);
				end

				a_cond   = zeros(nobs,m); 
				a        = zeros(nobs,m);
				P_cond   = zeros(nobs,m*m); 
				P        = zeros(nobs,m*m);
				a_smooth = zeros(nobs,m); 
				P_smooth = zeros(nobs,m*m);

				if timevar ~= 1 % then matrices are time invariant
				  Zt=Z;
				  dt=d;
				  Ht=H;
				  Tt=T;
				  ct=c;
				  %Rt=R;
				  Qt=Q;
				end

				for  i=1:nobs

				    yt=y(i,:)';                         % y(t) 

				    if timevar == 1
				      Zt=reshape(Z(i,:),n,m);           % Z(t) 
				      dt=d(i,:)';                       % d(t) 
				      Ht=reshape(H(i,:),n,n);           % H(t) 
				      Tt=reshape(T(i,:),m,m);           % T(t) 
				      ct=c(i,:)';                       % c(t) 
				      Rt=reshape(R(i,:),m,g);           % R(t) 
				      Qt=reshape(Q(i,:),g,g);           % Q(t) 
				    end;

				    % Prediction Equations 

				    at_1 = Tt*at + ct ;                  % a(t|t-1) formula(3.2.2a) 
				    Pt_1 = Tt*Pt*Tt' + Qt;%Rt*Qt*Rt' ;       % P(t|t-1) formula(3.2.2b) 

				    % Innovations 

				    yt_1 = Zt*at_1 + dt ;                % y(t|t-1) formula(3.2.18) 
				    vt = yt-yt_1 ;                       % v(t)     formula(3.2.19) 

				    % Updating Equations 

				    Ft = Zt*Pt_1*Zt' + Ht ;              % F(t)     formula(3.2.3c) 
				    inv_Ft = Ft\eye(size(Ft,1));         % Inversion de Ft                            
				    
				    at = at_1 + Pt_1*Zt'*inv_Ft*vt ;      % a(t)     formula(3.2.3a)   
				    Pt = Pt_1 - Pt_1*Zt'*inv_Ft*Zt*Pt_1 ; % P(t)     formula(3.2.3b)   

				    % Save results 

				    a(i,:)=at';
				    a_cond(i,:)=at_1';
				    P(i,:)=vecr(Pt)';
				    P_cond(i,:)=vecr(Pt_1)';

				end;

				% Smoothing Equations 

				a_smooth(nobs,:)=at';
				P_smooth(nobs,:)=vecr(Pt)';

				for i=(nobs-1):-1:1;
				    
				    if timevar == 1;
				      Tt=reshape(T(i+1,:),m,m);           % T(t+1) 
				    end;

				    Pt=reshape(P(i,:),m,m);               % P(t)     
				    Pt_1=reshape(P_cond(i+1,:),m,m);      % P(t+1|t) 
				    Pt_1_T=reshape(P_smooth(i+1,:),m,m);  % P(t+1|T) 
				    at=a(i,:)';                           % a(t)     
				    at_1=a_cond(i+1,:)';                  % a(t+1|t) 
				    at_1_T=a_smooth(i+1,:)';              % a(t+1|T) 

				    inv_Pt_1 = Pt_1\eye(size(Pt_1,1));         % Inversion de Ft           

				    P_star = Pt*Tt'*inv_Pt_1;

				    as = at + P_star*(at_1_T-at_1) ;          % a(t|T) formula(3.6.16a) 
				    Ps = Pt + P_star*(Pt_1_T-Pt_1)*P_star' ;  % P(t|T) formula(3.6.16b) 

				    %Added equation P(t|t-1)'
				    P2_smooth(i,:)=vecr((Pt_1_T*P_star')');
				    a_smooth(i,:)=as';
				    P_smooth(i,:)=vecr(Ps)';
				    
				end
			*/

            int nobs = y.R;
            int n = y.C;
            Matrix at = a0;
            Matrix Pt = P0;
            int m;
            int g;


            if (timeVar == 1)
            {
                m = Z.C / n;
                g = Q.C / m;
            }
            else
            {
                m = Z.C;
                g = Q.C;
            }

            Matrix a_cond = new Matrix(nobs, m);
            Matrix a = new Matrix(nobs, m);
            Matrix P_cond = new Matrix(nobs, m * m);
            Matrix P = new Matrix(nobs, m * m);
            Matrix a_smooth = new Matrix(nobs, m);
            Matrix P_smooth = new Matrix(nobs, m * m);
            Matrix P2_smooth = new Matrix(nobs - 1, 4);

            Matrix Zt = null;
            Matrix dt = null;
            Matrix Ht = null;
            Matrix Tt = null;
            Matrix ct = null;
            Matrix Qt = null;

            Matrix at_1 = null;
            Matrix Pt_1 = null;
            Matrix yt_1 = null;
            Matrix vt = null;
            Matrix inv_Ft = null;


            //% then matrices are time invariant
            if (timeVar != 1)
            {
                Zt = Z.Clone();
                dt = d.Clone();
                Ht = H.Clone();
                Tt = T.Clone();
                ct = c.Clone();
                Qt = Q.Clone();
            }

            for (int i = 0; i < nobs; i++)
            {
                Vector yt = y.GetRowReference(i);

                if (timeVar == 1)
                {
                    Zt = Z[i, Range.All];
                    Zt = Zt.Reshape(n, m);
                    dt = d[i, Range.All];
                    Ht = H[i, Range.All];
                    Ht = Ht.Reshape(n, n);

                    Tt = T[i, Range.All];
                    Tt = Tt.Reshape(m, m);
                    ct = c[i, Range.All].Transpose();
                    Qt = Q[i, Range.All];
                    Qt = Qt.Reshape(g, g);
                }

                //% Prediction Equations

                at_1 = Tt * at + ct;                  //% a(t|t-1) formula(3.2.2a)
                Pt_1 = Tt * Pt * Tt.Transpose() + Qt;    //%Rt*Qt*Rt' ;   //% P(t|t-1) formula(3.2.2b)

                //% Innovations

                yt_1 = Zt * at_1 + dt;                //% y(t|t-1) formula(3.2.18)
                vt = yt - yt_1;                       //% v(t)     formula(3.2.19)

                //% Updating Equations

                Matrix Ft = Zt * Pt_1 * Zt.Transpose() + Ht;              //% F(t)     formula(3.2.3c)
                inv_Ft = Ft.Inverse(); //Ft\Matrix.Identity(Ft.R);         //% Inversion de Ft

                at = at_1 + Pt_1 * Zt.Transpose() * inv_Ft * vt;      //% a(t)     formula(3.2.3a)
                Pt = Pt_1 - Pt_1 * Zt.Transpose() * inv_Ft * Zt * Pt_1; //% P(t)     formula(3.2.3b)

                //% Save results
                a[i, Range.All] = at.Transpose();
                a_cond[i, Range.All] = at_1.Transpose();
                P[i, Range.All] = Vecr(Pt).Transpose();
                P_cond[i, Range.All] = Vecr(Pt_1).Transpose();
            }

            /*
            //% Smoothing Equations
            */
            a_smooth[nobs - 1, Range.All] = at.Transpose();
            P_smooth[nobs - 1, Range.All] = Vecr(Pt).Transpose();

            for (int i = (nobs - 2); i >= 0; i--)
            {
                if (timeVar == 1)
                {
                    Tt = T[i + 1, Range.All].Reshape(m, m);
                }

                Pt = P[i, Range.All].Reshape(m, m);

                Pt_1 = P_cond[i + 1, Range.All].Reshape(m, m);

                Matrix Pt_1_T = P_smooth[i + 1, Range.All];
                Pt_1_T = Pt_1_T.Reshape(m, m);
                at = a[i, Range.All].Transpose();
                at_1 = a_cond[i + 1, Range.All].Transpose();
                Matrix at_1_t = a_smooth[i + 1, Range.All].Transpose();
                Matrix inv_Pt_1 = Pt_1.Inverse();
                Matrix P_star = Pt * Tt.Transpose() * inv_Pt_1;
                Matrix a_s = at + P_star * (at_1_t - at_1);
                Matrix Ps = Pt + P_star * (Pt_1_T - Pt_1) * P_star;

                P2_smooth[i, Range.All] = Vecr((Pt_1_T * P_star.Transpose()).Transpose()).Transpose();
                a_smooth[i, Range.All] = a_s.Transpose();
                P_smooth[i, Range.All] = Vecr(Ps).Transpose();
            }

            return new Matrix[] { a_smooth, P_smooth, P2_smooth };
        }

        public static Matrix[] KalmanFilter(Matrix y, Matrix Z, Matrix d, Matrix T, Matrix c, Matrix a0,
                             Matrix P0, Matrix H, Matrix Q, int timeVar)
        {
            /*
                function [a_smooth,P_smooth,P2_smooth] = Kalman_Smoothing(y,Z,d,T,c,a0,P0,H,Q,timevar)
                % Kalman smoothing
                % code is a Matlab translation of Thierry Roncalli's.

                nobs = size(y,1); 
                n    = size(y,2); 
                at   = a0; 
                Pt   = P0; 


                if timevar == 1
                    m=size(Z,2)/n; 
                    g=size(Q,2)/m;
                else
                    m=size(Z,2); 
                    g=size(Q,2);
                end

                a_cond   = zeros(nobs,m); 
                a        = zeros(nobs,m);
                P_cond   = zeros(nobs,m*m); 
                P        = zeros(nobs,m*m);
                a_smooth = zeros(nobs,m); 
                P_smooth = zeros(nobs,m*m);

                if timevar ~= 1 % then matrices are time invariant
                  Zt=Z;
                  dt=d;
                  Ht=H;
                  Tt=T;
                  ct=c;
                  %Rt=R;
                  Qt=Q;
                end

                for  i=1:nobs

                    yt=y(i,:)';                         % y(t) 

                    if timevar == 1
                      Zt=reshape(Z(i,:),n,m);           % Z(t) 
                      dt=d(i,:)';                       % d(t) 
                      Ht=reshape(H(i,:),n,n);           % H(t) 
                      Tt=reshape(T(i,:),m,m);           % T(t) 
                      ct=c(i,:)';                       % c(t) 
                      Rt=reshape(R(i,:),m,g);           % R(t) 
                      Qt=reshape(Q(i,:),g,g);           % Q(t) 
                    end;

                    % Prediction Equations 

                    at_1 = Tt*at + ct ;                  % a(t|t-1) formula(3.2.2a) 
                    Pt_1 = Tt*Pt*Tt' + Qt;%Rt*Qt*Rt' ;       % P(t|t-1) formula(3.2.2b) 

                    % Innovations 

                    yt_1 = Zt*at_1 + dt ;                % y(t|t-1) formula(3.2.18) 
                    vt = yt-yt_1 ;                       % v(t)     formula(3.2.19) 

                    % Updating Equations 

                    Ft = Zt*Pt_1*Zt' + Ht ;              % F(t)     formula(3.2.3c) 
                    inv_Ft = Ft\eye(size(Ft,1));         % Inversion de Ft                            
				    
                    at = at_1 + Pt_1*Zt'*inv_Ft*vt ;      % a(t)     formula(3.2.3a)   
                    Pt = Pt_1 - Pt_1*Zt'*inv_Ft*Zt*Pt_1 ; % P(t)     formula(3.2.3b)   

                    % Save results 

                    a(i,:)=at';
                    a_cond(i,:)=at_1';
                    P(i,:)=vecr(Pt)';
                    P_cond(i,:)=vecr(Pt_1)';

                end;

                % Smoothing Equations 

                a_smooth(nobs,:)=at';
                P_smooth(nobs,:)=vecr(Pt)';

                for i=(nobs-1):-1:1;
				    
                    if timevar == 1;
                      Tt=reshape(T(i+1,:),m,m);           % T(t+1) 
                    end;

                    Pt=reshape(P(i,:),m,m);               % P(t)     
                    Pt_1=reshape(P_cond(i+1,:),m,m);      % P(t+1|t) 
                    Pt_1_T=reshape(P_smooth(i+1,:),m,m);  % P(t+1|T) 
                    at=a(i,:)';                           % a(t)     
                    at_1=a_cond(i+1,:)';                  % a(t+1|t) 
                    at_1_T=a_smooth(i+1,:)';              % a(t+1|T) 

                    inv_Pt_1 = Pt_1\eye(size(Pt_1,1));         % Inversion de Ft           

                    P_star = Pt*Tt'*inv_Pt_1;

                    as = at + P_star*(at_1_T-at_1) ;          % a(t|T) formula(3.6.16a) 
                    Ps = Pt + P_star*(Pt_1_T-Pt_1)*P_star' ;  % P(t|T) formula(3.6.16b) 

                    %Added equation P(t|t-1)'
                    P2_smooth(i,:)=vecr((Pt_1_T*P_star')');
                    a_smooth(i,:)=as';
                    P_smooth(i,:)=vecr(Ps)';
				    
                end
            */

            int nobs = y.R;
            int n = y.C;
            int m = 0;
            int g = 0;
            Matrix at = a0;
            Matrix Pt = P0;
            Matrix logl = new Matrix(nobs, 1);

            if (timeVar == 1)
            {
                m = Z.C / n;
                g = Q.C / n;
            }
            else
            {
                m = Z.C;
                g = Q.C;
            }

            Matrix y_cond = new Matrix(nobs, n);
            Matrix v = new Matrix(nobs, n);
            Matrix a_cond = new Matrix(nobs, m);
            Matrix a = new Matrix(nobs, m);
            Matrix P_cond = new Matrix(nobs, m * m);
            Matrix P = new Matrix(nobs, m * m);
            Matrix F = new Matrix(nobs, n * n);

            Matrix Zt = new Matrix();
            Matrix dt = new Matrix();
            Matrix Ht = new Matrix();
            Matrix Tt = new Matrix();
            Matrix ct = new Matrix();
            Matrix Qt = new Matrix();

            if (timeVar != 1)
            {
                Zt = Z;
                dt = d;
                Ht = H;
                Tt = T;
                ct = c;
                Qt = Q;
            }

            for (int i = 0; i < nobs; i++)
            {
                Matrix yt = y[i, Range.All];
                yt = yt.Transpose();

                if (timeVar == 1)
                {
                    Zt = Z[i, Range.All];
                    Zt = Zt.Reshape(n, m);
                    dt = d[i, Range.All];
                    dt = dt.Transpose();
                    Ht = H[i, Range.All];
                    Ht = Ht.Reshape(n, n);
                    Tt = T[i, Range.All];
                    Tt = Tt.Reshape(m, m);
                    ct = c[i, Range.All];
                    ct = ct.Transpose();
                    Qt = Q[i, Range.All];
                    Qt = Qt.Reshape(g, g);
                }

                //% Prediction Equations
                Matrix at_1 = Tt * at + ct;
                Matrix Pt_1 = Tt * Pt * Tt.Transpose() + Qt;

                //% Innovations
                Matrix yt_1 = Zt * at_1 + dt;
                Matrix vt = yt - yt_1;

                //% Updating Equation
                Matrix Ft = Zt * Pt_1 * Zt.Transpose() + Ht;
                Matrix inv_Ft = Ft.Inverse();

                at = at_1 + Pt_1 * Zt.Transpose() * inv_Ft * vt;
                Pt = Pt_1 - Pt_1 * Zt.Transpose() * inv_Ft * Zt * Pt_1;

                //% Save results
                y_cond[i, Range.All] = yt_1.Transpose();
                v[i, Range.All] = vt.Transpose();
                a[i, Range.All] = at.Transpose();
                a_cond[i, Range.All] = at_1.Transpose();
                P[i, Range.All] = Vecr(Pt).Transpose();
                P_cond[i, Range.All] = Vecr(Pt_1).Transpose();
                F[i, Range.All] = Vecr(Ft).Transpose();

                double dFt = Ft.Determinant();
                if (dFt <= 0)
                    dFt = 1e-10;


                // The matrix here is 1x1.
                logl[i, 0] = (double)(-(n / 2) * Math.Log(2 * Math.PI) - 0.5 * Math.Log(dFt) - 0.5 * vt.Transpose() * inv_Ft * vt);
            }

            return new Matrix[] { y_cond, v, a, a_cond, P, P_cond, F, logl };
        }

        public static Matrix[] FuOptimization(Matrix y, Matrix Hhat, Matrix dhat, Matrix Ahat, Matrix Chat, Matrix x0, Matrix W0, Matrix Qhat,
                                          Matrix What, int tol, int T, double dt, Vector maturities, int nf)
        {
            //function [Hhat1,dhat1,Ahat1,Chat1,Qhat1,What1,ll1,estimate]=fuopt(y,Hhat,...
            //dhat,Ahat,Chat,x0,W0,Qhat,What,tol,T,dt,matur,nf)

            int iter = 0;
            Matrix ll = new Matrix(new double[,] { { 1, 1 + 2 * tol } });
            Matrix estimate = new Matrix(1, 7);
            Matrix a_smooth = new Matrix();
            Matrix P_smooth = new Matrix();
            Matrix P2_smooth = new Matrix();

            while (Math.Abs(ll[0, iter + 1] - ll[0, iter]) > tol)
            {
                /*[a_smooth,P_smooth,P2_smooth] = Kalman_Smoothing(y,Hhat,dhat,Ahat,Chat,x0,W0,Qhat,What,0);*/
                Matrix[] kalOutput = KalmanSmoothing(y, Hhat, dhat, Ahat, Chat, x0, W0, Qhat, What, 0);
                a_smooth = kalOutput[0];
                P_smooth = kalOutput[1];
                P2_smooth = kalOutput[2];

                Matrix mu = a_smooth.Transpose();

                Matrix sumP = P_smooth.ColumnSum().Reshape(2, 2).Transpose();
                Matrix sumP1 = P_smooth[new Range(0, P_smooth.R - 2), Range.All].ColumnSum().Reshape(2, 2).Transpose();
                Matrix sumP2 = P_smooth[new Range(1, P_smooth.R - 1), Range.All].ColumnSum().Reshape(2, 2).Transpose();
                Matrix sumP21 = P2_smooth.ColumnSum().Reshape(2, 2).Transpose();

                dhat = DM(y.Transpose(), mu, sumP);
                Hhat = HM(y.Transpose(), mu, sumP);
                Qhat = QM(y.Transpose(), dhat, Hhat, mu, sumP, T);
                Chat = CM(mu, sumP21, sumP1);
                Ahat = AM(mu, sumP21, sumP1);
                What = WM(mu, sumP21, sumP2, sumP1, Ahat, Chat, T);

                /*
	                //%Calculate mu, mu[1], and mu[2], summu,summu1,summu2
	                mu=a_smooth';

	                //%Calculate SumP, SumP[1], SumP[2],SumP[2,1]
	                sumP=reshape(sum(P_smooth),2,2)';
	                sumP1=reshape(sum(P_smooth(1:end-1,:)),2,2)';
	                sumP2=reshape(sum(P_smooth(2:end,:)),2,2)';
	                sumP21=reshape(sum(P2_smooth),2,2)';

	                //%Update Matrices
	                dhat=dm(y.Traspose(),mu,sumP);
	                Hhat=hm(y.Traspose(),mu,sumP);
	                Qhat=qm(y.Traspose(),dhat,Hhat,mu,sumP,T);
	                Chat=cm(mu,sumP21,sumP1);
	                Ahat=am(mu,sumP21,sumP1);
	                What=wm(mu,sumP21,sumP2,sumP1,Ahat,Chat,T);
                */


                /*
	                //%Compute the likelihood of Y under estimated parameters
	                [y_cond,v,a,a_cond,P,P_cond,F,logl]=kalman_filter(y,Hhat,dhat,Ahat,Chat,x0,W0,Qhat,What,0);
	                estimate(iter,:)=extract(Hhat,dhat,Ahat,Chat,What,dt,matur,nf);
	                theta2=constraints(estimate(iter,:));
                */
                kalOutput = KalmanFilter(y, Hhat, dhat, Ahat, Chat, x0, W0, Qhat, What, 0);
                Matrix y_cond = kalOutput[0];
                Matrix v = kalOutput[1];
                Matrix a = kalOutput[2];
                Matrix a_cond = kalOutput[3];
                Matrix P = kalOutput[4];
                Matrix P_cond = kalOutput[5];
                Matrix F = kalOutput[6];
                Matrix logl = kalOutput[7];

                estimate[iter, Range.All] = Extract(Hhat, dhat, Ahat, Chat, What, dt, maturities, nf);
                Console.WriteLine(estimate.ToString());
                Vector theta2 = Constraints((Vector)estimate[iter, Range.All]);
                Console.WriteLine(theta2.ToString());

                Matrix[] outputMatrices = Matrices(theta2[0], theta2[1], theta2[2],
                                                    theta2[3], theta2[4], theta2[5],
                                                    theta2[6], Qhat, maturities, dt);
                Chat = outputMatrices[0];
                Ahat = outputMatrices[1];
                What = outputMatrices[2];
                dhat = outputMatrices[3];
                Hhat = outputMatrices[4];
                Qhat = outputMatrices[5];

                iter = iter + 1;
                // TODO is this the first element in practice?
               
                ll[0, iter + 1] = logl[Range.All, 0].Sum();
            }


            return new Matrix[] { Hhat, dhat, Ahat, Chat, Qhat, What, ll, estimate };

        }

        public static Matrix Extract(Matrix Hhat, Matrix dhat, Matrix Ahat, Matrix Chat, Matrix What, double dt, Vector maturities, int nf)
        {
            /*
				function y=extract(Hhat,dhat,Ahat,Chat,What,dt,matur,nf)

				%k=-(log(Ahat(1,1))+sum(log(Hhat(:,1))))/(dt+sum(matur));
				%mue=Chat(2,1)/dt;
				%sigmax=sqrt(What(1,1)/(1-exp(-2*k*dt))*(2*k));
				%sigmae=sqrt(What(2,2)/dt);
				%pxe=(What(1,2)+What(2,1))*k/(((1-exp(-k*dt))*sigmax*sigmae)*2);

				k=-(log(abs(Ahat(1,1)))+sum(log(abs(Hhat(:,1)))))/(dt+sum(matur));
				mue=Chat(2,1)/dt;
				sigmax=sqrt(abs(What(1,1))/(1-exp(-2*k*dt))*(2*k));
				sigmae=sqrt(abs(What(2,2))/dt);
				pxe=(abs(What(1,2))+abs(What(2,1)))*k/(((1-exp(-k*dt))*sigmax*sigmae)*2);

				for i=1:nf
				p1=(1-exp(-2*k*matur(i)))*(sigmax)^2/(2*k);
				p2=(sigmae)^2*matur(i);
				p3=2*(1-exp(-k*matur(i)))*pxe*sigmax*sigmae/k;
				dd(i,1)=dhat(i,1)-.5*(p1+p2+p3)-mue*matur(i);

				aa(i,1)=-matur(i);
				aa(i,2)=-(1-exp(-k*matur(i)))/k;
				end

				value=inv((aa'*aa))*aa'*dd;
				lambdae=value(1);
				lambdax=value(2);

				y=[k,mue,sigmax,sigmae,pxe,lambdae,lambdax];
            */

            double k = -((Matrix.Log(Matrix.Abs(new Matrix(new[] { Ahat[0, 0] }))) + (Matrix.Log(Matrix.Abs(Hhat[Range.All, 0]))).ColumnSum()) / (dt + maturities.Sum()))[0, 0];
            double mue = Chat[1, 0] / dt;
            double sigmax = Math.Sqrt(Math.Abs(What[0, 0]) / (1 - Math.Exp(-2 * k * dt)) * (2 * k));
            double sigmae = Math.Sqrt(Math.Abs(What[1, 1]) / dt);
            double pxe = (Math.Abs(What[0, 1]) + Math.Abs(What[1, 0])) * k / (((1 - Math.Exp(-k * dt)) * sigmax * sigmae) * 2);
            Matrix dd = new Matrix(nf, 1);
            Matrix aa = new Matrix(nf, 2);
            for (int i = 0; i < nf; i++)
            {
                double p1 = (1 - Math.Exp(-2 * k * maturities[i])) * Math.Pow(sigmax, 2) / (2 * k);
                double p2 = Math.Pow(sigmae, 2) * maturities[i];
                double p3 = 2 * (1 - Math.Exp(-k * maturities[i])) * pxe * sigmax * sigmae / k;
                dd[i, 0] = dhat[i, 0] - 0.5 * (p1 + p2 + p3) - mue * maturities[i];
                aa[i, 0] = -maturities[i];
                aa[i, 1] = -(1 - Math.Exp(-k * maturities[i])) / k;
            }

            Matrix values = (aa.Transpose() * aa).Inverse() * aa.Transpose() * dd;
            double lambdae = values[0, 0];
            double lambdax = values[1, 0];

            return new Matrix(new double[,] { { k, mue, sigmax, sigmae, pxe, lambdae, lambdax } });
        }

        public static Vector Constraints(Vector x)
        {
            /*
				function y=constraints(x) %theta=[.8,.05,.2,.15,.5,.03,.15];
				%[k,mue,sigmax,sigmae,pxe,lambdae,lambdax]
				%data theta0=[1.49,.03,.286,.145,.3,.0415,.157];
				%simulation theta=[1.2,.1,.4,.15,.2,.05,.05];
				%c=[.5,1;.01,.5;0,1;0,1;0,1;0,.1;0,.1];
				c=[1,2;.01,.2;0,.5;0,.5;.1,.4;0,.2;0,.2];
				[m,n]=size(c);

				for i=1:m
				    
				if (x(i)<c(i,1) || x(i)>c(i,2))
				x(i)=c(i,1) + rand*(c(i,2)-c(i,1));
				else
				x(i)=x(i);
				end

				end
				y=x;
			*/

            Vector y = x.Clone();

            Matrix c = new Matrix(new double[,] { { 1, 2 }, { .01, .2 }, { 0, .5 }, { 0, .5 }, { .1, .4 }, { 0, .2 }, { 0, .2 } });
            int m = c.R;

            Random r = new Random();
            for (int i = 0; i < m; i++)
            {
                if (y[i] < c[i, 0] || y[i] > c[i, 1])
                {

                    y[i] = c[i, 0] + r.NextDouble() * (c[i, 1] - c[i, 0]);
                }
                else
                {
                    // TODO: why?
                    y[i] = y[i];
                }
            }

            return y;
        }

        public static Matrix FutureErrors(Vector maturities, Matrix y, Matrix dhat, Matrix Hhat, Matrix a1)
        {
            /*
                function fe = FuturesError(matur,y,dhat,Hhat,a1)
                for i=1:length(matur)
                est=dhat(i)+Hhat(i,1)*a1(:,1)+Hhat(i,2)*a1(:,2);
                obs=y(i,:)';
                fe(i,2)=norm(est-obs,'fro')/norm(obs,'fro');
                end
			 */

            Matrix fe = new Matrix();
            for (int i = 0; i < maturities.Count; i++)
            {
                //Matrix est = dhat[i]+Hhat[i,1]*a1[Range.All, 1]+Hhat[i,2]*a1[Range.All,2];
                Matrix obs = y[i, Range.All];
                obs = obs.Transpose();
                //fe[i,2] =  TODO NORM!!!
            }

            return fe;
        }

        public static Matrix[] Matrices(double k, double mue, double sigmax, double sigmae, double pxe, double lambdae, double lambdax, Matrix s, Vector maturities, double dt)
        {
            /* function [C,A,W,d,H,Q]=matrices(k,mue,sigmax,sigmae,pxe,lambdae,lambdax,s,matur,dt)
                %maturity times of futures contracts
                nf=length(matur);

                %Value of C
                C=[0;mue*dt];

                % Value of A
                A=[exp(-k*dt),0;0,1];

                %Value of W
                xx=(1-exp(-2*k*dt))*(sigmax)^2/(2*k);
                xy=(1-exp(-k*dt))*pxe*sigmax*sigmae/k;
                yx=(1-exp(-k*dt))*pxe*sigmax*sigmae/k;
                yy=(sigmae)^2*dt;
                W=[xx,xy;yx,yy];

                %Value of d and H
                for z=1:nf
                %p1, p2, and p3 are components of A(T) Equation 9
                p1=(1-exp(-2*k*matur(z)))*(sigmax)^2/(2*k);
                p2=(sigmae)^2*matur(z);
                p3=2*(1-exp(-k*matur(z)))*pxe*sigmax*sigmae/k;
                d(z,1)=(mue-lambdae)*matur(z)-(1-exp(-k*matur(z)))*lambdax/k+.5*(p1+p2+p3); %Value of d(t)
                H(z,1)=exp(-k*matur(z)); % Value of F column 1
                H(z,2)=1;  % Value of F column 2
                end

                [m,n]=size(s);
                if n==1
                %Value of Q
                Q=diag(s);
                else
                Q=s;
                end
            */

            //%maturity times of futures contracts
            int nf = maturities.Count;

            //%Value of C
            Vector C = new Vector { 0, mue * dt };

            //% Value of A
            Matrix A = new Matrix(new[,] { { Math.Exp(-k * dt), 0 }, { 0, 1 } });

            //%Value of W
            double xx = (1 - Math.Exp(-2 * k * dt)) * Math.Pow(sigmax, 2) / (2 * k);
            double xy = (1 - Math.Exp(-k * dt)) * pxe * sigmax * sigmae / k;
            double yx = (1 - Math.Exp(-k * dt)) * pxe * sigmax * sigmae / k;
            double yy = Math.Pow(sigmae, 2) * dt;
            Matrix W = new Matrix(new[,] { { xx, xy }, { yx, yy } });
            Vector d = new Vector(nf);
            Matrix H = new Matrix(nf, 2);

            //%Value of d and H
            for (int z = 0; z < nf; z++)
            {
                //%p1, p2, and p3 are components of A(T) Equation 9
                double p1 = (1 - Math.Exp(-2 * k * maturities[z])) * Math.Pow(sigmax, 2) / (2 * k);
                double p2 = Math.Pow(sigmae, 2) * maturities[z];
                double p3 = 2 * (1 - Math.Exp(-k * maturities[z])) * pxe * sigmax * sigmae / k;
                //%Value of d(t)
                d[z] = (mue - lambdae) * maturities[z] - (1 - Math.Exp(-k * maturities[z])) * lambdax / k + .5 * (p1 + p2 + p3);
                //% Value of F column 1
                H[z, 0] = Math.Exp(-k * maturities[z]);
                //% Value of F column 2
                H[z, 1] = 1;
            }

            int n = s.C;
            Matrix Q;
            if (n == 1)
            {
                //%Value of Q
                Q = Matrix.Diag(s);
            }
            else
            {
                Q = s;
            }

            return new Matrix[] { C, A, W, d, H, Q };
        }

        // %This program implement the crude-oil dataset used in the schwartz-smith paper.
        public static void Main(string[] args)
        {
            Matrix x = new Matrix(DVPLI.ArrayHelper.FromTextFile("crudeoilspot.txt"));
            Matrix y = new Matrix(DVPLI.ArrayHelper.FromTextFile("crudeoil.txt")).Transpose();

            Vector maturities = new Vector { 1.0 / 12, 5.0 / 12, 9.0 / 12, 13.0 / 12, 17.0 / 12 };
            int nf = y.R;
            int T = y.C;

            double dt = 7.0 / 360;

            Matrix x0 = new Matrix(new double[,] { { .119 }, { 2.857 } }); //Initial state vector m(t)=E[xt;et]
            Matrix W0 = new Matrix(new double[,] { { .1, .1 }, { .1, .1 } }); //%Inital covariance matrix C(t)=cov[xt,et]

            //%% Test Smith&Schwartz original TODO: ADD TO TEST CASE LATER?
            //%Tolerence level for the EM Estimation of dataset
            double l = TestOriginalSS(dt, y, maturities, x.Transpose(), x0, W0);
            int tol = 1;

            //%% Initial Guess at Parameters for em estimation
            //%(k,mue,sigmax,sigmae,pxe,lambdae,lambdax)
            //%theta=[1,.05,.1,.1,.5,0,.1];

            Vector Theta = new Vector();
            Random random = new Random();
            for (int i = 0; i < 7; i++)
            {
                Theta.Add(Math.Abs(random.NextDouble()));
            }
            Theta = new Vector { 1, .05, .1, .1, .5, 0, .1 };

            Vector s2 = new Vector { .005, .005, .005, .005, .005 };

            //%% Generate the matrices that correspond the given parameters
            Matrix[] outputMat = Matrices(Theta[0], Theta[1], Theta[2], Theta[3], Theta[4], Theta[5], Theta[6], s2, maturities, dt);
            Matrix Chat = outputMat[0];
            Matrix Ahat = outputMat[1];
            Matrix What = outputMat[2];
            Matrix dhat = outputMat[3];
            Matrix Hhat = outputMat[4];
            Matrix Qhat = outputMat[5];

            //%% Fully unstructured optimization
            outputMat = FuOptimization(y.Transpose(), Hhat, dhat, Ahat, Chat, x0, W0, Qhat, What, tol, T, dt, maturities, nf);

            Matrix Hhat1 = outputMat[0];
            Matrix dhat1 = outputMat[1];
            Matrix Ahat1 = outputMat[2];
            Matrix Chat1 = outputMat[3];
            Matrix Qhat1 = outputMat[4];
            Matrix What1 = outputMat[5];
            Matrix ll1 = outputMat[6];
            Matrix e1 = outputMat[7];

            /*
            //%% Fully unstructured optimization
            [Hhat1,dhat1,Ahat1,Chat1,Qhat1,What1,ll1,e1]=fuopt(y',Hhat,dhat,Ahat,Chat,...
                x0,W0,Qhat,What,tol,T,dt,matur,nf);

            //%% Estimate the hidden state given the estimates of the parameters
            [y_cond1,v1,a1,a_cond1,P1,P_cond1,F1,logl1]=kalman_filter(y',Hhat1,dhat1,...
                Ahat1,Chat1,x0,W0,Qhat1,What1,0);

            //%% Plot the estimated data and observed data given by the Fully unstructured
            //%% optimization
            GraphComparisonResult(a1,x,ll1)

            thetahat(1,:)=extract(Hhat1,dhat1,Ahat1,Chat1,What1,dt,matur,nf);
            */

            outputMat = KalmanFilter(y.Transpose(), Hhat1, dhat1, Ahat1, Chat1, x0, W0, Qhat1, What1, 0);
            Matrix y_cond1 = outputMat[0];
            Matrix v1 = outputMat[1];
            Matrix a1 = outputMat[2];
            Matrix a_cond1 = outputMat[3];
            Matrix P1 = outputMat[4];
            Matrix P_cond1 = outputMat[5];
            Matrix F1 = outputMat[6];
            Matrix logl1 = outputMat[7];

            Matrix thetahat = new Matrix(1, 7);
            thetahat[0, Range.All] = Extract(Hhat1, dhat1, Ahat1, Chat1, What1, dt, maturities, nf);

            //%% Futures price errors
            Matrix fe = FutureErrors(maturities, y, dhat1, Hhat1, a1);
            Console.WriteLine(thetahat);
            //thetahat
        }
    }

    internal static class MatrixExtensions
    {
        internal static Matrix ColumnSum(this Matrix t)
        {
            Matrix output = new Matrix(1, t.C, false);
            for (int c = 0; c < t.C; c++)
            {
                double sum = 0;
                for (int r = 0; r < t.R; r++)
                    sum += t[r, c];

                output[0, c] = sum;
            }

            return output;
        }
    }
}

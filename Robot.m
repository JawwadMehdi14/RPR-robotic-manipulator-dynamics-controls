classdef Robot
    properties
        DH_table;
        TransoformationMatrix;
        Joints;
        GeometricalJacobianMatrix;
        AnalyticalJacobianMatrix;
        AnalyticalJacobianMatrixComplete;
    end
    properties 
        r1,inertia_vars,q, q_dot, q_dot_dot, masses, lengths, g, f_e, mu_e, P_COM, phi, x_dot_dot, Fs, Fv;
    end
    methods
        %%Default constructor
        function obj = Robot()
            obj.r1 = sym("r1",'real'); %base frame r in z axis
            obj.g = sym("g",'real');
            obj.q = [sym("q1",'real'); sym("q2",'real'); sym("q3",'real')];
            obj.q_dot = [sym("q1_dot",'real'); sym("q2_dot",'real'); sym("q3_dot",'real')];
            obj.q_dot_dot = [sym("q1_dot_dot",'real'); sym("q2_dot_dot",'real'); sym("q3_dot_dot",'real')];
            obj.masses = [sym("m1",'real'), sym("m2",'real'), sym("m3",'real')];
            obj.lengths = [[sym("a1",'real'),sym("b1",'real'),sym("c1",'real')];
                            [sym("a2",'real'),sym("b2",'real'),sym("c2",'real')];
                            [sym("a3",'real'),sym("b3",'real'),sym("c3",'real')]];
            obj.phi = [ sym("phi1",'real'), sym("phi2",'real'), sym("phi3",'real')];
            obj.x_dot_dot = [sym("x_dot_dot",'real');sym("y_dot_dot",'real'); sym("z_dot_dot",'real'); sym("phi1_dot_dot",'real');sym("phi2_dot_dot",'real');sym("phi3_dot_dot",'real') ];
            obj.P_COM = [
                [-obj.lengths(1,3)/2,0,0];
                [0,0,-obj.lengths(2,3)/2];
                [0,0,-obj.lengths(3,3)/2];];
            obj.Joints = ["Revolut","Prismatic","Prismatic"];
            obj.DH_table = [
                % ["theta", "alpha","a","d"]
                [0,pi/2, 0, obj.r1];
                [obj.q(1),0,obj.lengths(1,3),0];
                [pi/2,pi/2,0,obj.lengths(2,3)+obj.q(2)];
                [-pi/2,0,0,obj.lengths(3,3)+obj.q(3)];
                ];
             obj.f_e = [sym("fex",'real');sym("fey",'real');sym("fez",'real')];
             obj.mu_e = [sym("muex",'real');sym("muey",'real');sym("muez",'real')];
             obj.Fs = [sym("Fs1","real"),sym("Fs2","real"),sym("Fs3","real")];
             obj.Fv = [sym("Fv1","real"),sym("Fv2","real"),sym("Fv3","real")];
             obj.TransoformationMatrix = obj.getTransformationMatrix(obj.DH_table, 4);
             obj.GeometricalJacobianMatrix = obj.generateJacobianMatrix(obj.DH_table,obj.Joints, 3);
             obj.AnalyticalJacobianMatrix = obj.generateAnalyticalJacobianMatrix(obj.TransoformationMatrix(1:3,4),obj.q);
             [~,obj.AnalyticalJacobianMatrixComplete] = obj.generateAnalyticalJacobianComplete(obj.GeometricalJacobianMatrix);
        end

        function kinematics = directKinematics(obj, t_1, d_1, d_2)
           kinematics = double(subs(obj.TransoformationMatrix, [obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3),obj.q(1), obj.q(3),obj.q(2)], [0.4, 0.4, 0.46, 0.15, t_1,d_2, d_1]));
        end

        function [q1_ik_val, q2_ik_val, q3_ik_val] = inverseKinematics(obj, t_1, d_1, d_2)

            syms q1 q2 q3 real

            a1 = obj.r1;
            a2 = obj.lengths(1,3);
            a3 = obj.lengths(2,3);
            a4 = obj.lengths(3,3);

            Px = cos(q1)*(a4 + q3) + a2*cos(q1);
            Py = -a3 - q2;
            Pz = a1 + sin(q1)*(a4 + q3) + a2*sin(q1);
        
            Px_val = (subs(Px, [q1, q2, q3], [t_1, d_1, d_2]));
            Py_val = (subs(Py, [q1, q2, q3], [t_1, d_1, d_2]));
            Pz_val = (subs(Pz, [q1, q2, q3], [t_1, d_1, d_2]));
        
            % Inverse kinematics
            q2_ik_val = -a3 - Py_val;
            q3_ik_val = (Px_val / cos(t_1)) - a4 - a2;
            q1_ik_val = acos(Px_val / (a4 + q3_ik_val + a2));
        end

        function JacobianMatrixSolved = solveGeometricalJacobianMatrix(obj, t_1,d_1,d_2)
            % JacobianMatrixSolved = double(subs(obj.GeometricalJacobianMatrix,[obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3), obj.q(1),obj.q(2),obj.q(3)],[0.4, 0.4, 0.46, 0.15,t_1,d_1,d_2]));
            JacobianMatrixSolved = double(subs(obj.GeometricalJacobianMatrix,[obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3), obj.q(1),obj.q(2),obj.q(3)],[0.4, 0.4, 0.46, 0.15,t_1,d_1,d_2]));
        end

        function AnalyticalJacobianMatrixSolved = solveAnalyticalJacobianMatrix(obj, t_1,d_1,d_2)
            %AnalyJac = obj.AnalyticalJacobianMatrix;
            AnalyticalJacobianMatrixSolved = double(subs(obj.AnalyticalJacobianMatrix,[obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3), obj.q(1),obj.q(2),obj.q(3)],[0.4, 0.4, 0.46, 0.15,t_1,d_1,d_2]));
        end

        function torque = newtonEulerEquationOfMotion(obj)
            [w, w_dot, p_dot_dot_c] = obj.forwardRecursion(obj.q, obj.q_dot, obj.q_dot_dot, obj.g, obj.Joints, obj.lengths);
            torque = obj.backwardEquation(obj.masses,p_dot_dot_c, w, w_dot, obj.f_e, obj.mu_e, obj.Joints, obj.lengths );
        end

        function test(obj, robot, q_values)

            r1_value = 0.4;
            density = 2700; %Almunium density
            mass_values = [pi * obj.lengths(1,1)^2*obj.lengths(1,3)*density, obj.lengths(2,1)*obj.lengths(2,2)*obj.lengths(2,3)*density, obj.lengths(3,1)*obj.lengths(3,2)*obj.lengths(3,3)*density];
            Fs_values = [1.5,1.5,1.5]*0;
            Fv_values = [0.5,0.5,0.5]*0;
            lenghts_values = [
               [0.02, 2*pi*0.02, 0.4]; 
               [0.03, 0.03, 0.46]; 
               [0.01, 0.01, 0.15]
               ];
            mass_values = subs(mass_values, [obj.lengths(1,1),obj.lengths(1,2), obj.lengths(1,3), obj.lengths(2,1),obj.lengths(2,2), obj.lengths(2,3),obj.lengths(3,1),obj.lengths(3,2), obj.lengths(3,3)], [lenghts_values(1,1), lenghts_values(1,2), lenghts_values(1,3),lenghts_values(2,1), lenghts_values(2,2), lenghts_values(2,3),lenghts_values(3,1), lenghts_values(3,2), lenghts_values(3,3)]);
            q_dot_values = [0.5,-0.3,0.1]*0;
            q_dot_dot_values = [0.2,0.1,-0.2]*0;
            g_value = 9.806;
            f_e_values = [0;0;0];
            mu_e_values = [0;0;0];
            
            parameters_old = [obj.masses(1), obj.masses(2), obj.masses(3),obj.g, obj.q(1), obj.q(2), obj.q(3),obj.r1, obj.lengths(1,1),obj.lengths(1,2), obj.lengths(1,3), obj.lengths(2,1),obj.lengths(2,2), obj.lengths(2,3),obj.lengths(3,1),obj.lengths(3,2), obj.lengths(3,3), obj.q_dot(1),obj.q_dot(2),obj.q_dot(3),obj.q_dot_dot(1),obj.q_dot_dot(2),obj.q_dot_dot(3), obj.f_e(1),obj.f_e(2),obj.f_e(3),obj.mu_e(1),obj.mu_e(2),obj.mu_e(3), obj.Fs, obj.Fv];
            parameters_new  = [mass_values(1), mass_values(2), mass_values(3), g_value, q_values(1),q_values(2),q_values(3),r1_value, lenghts_values(1,1), lenghts_values(1,2), lenghts_values(1,3),lenghts_values(2,1), lenghts_values(2,2), lenghts_values(2,3),lenghts_values(3,1), lenghts_values(3,2), lenghts_values(3,3),q_dot_values(1),q_dot_values(2),q_dot_values(3),q_dot_dot_values(1),q_dot_dot_values(2),q_dot_dot_values(3), f_e_values(1),f_e_values(2),f_e_values(3), mu_e_values(1),mu_e_values(2), mu_e_values(3), Fs_values, Fv_values];
            
            disp("q vector: ");
            disp(q_values);
            
            config = homeConfiguration(robot);
            config(1).JointPosition = q_values(1);
            config(2).JointPosition = q_values(2);
            config(3).JointPosition = q_values(3);
            
            transformResultToolKit = getTransform(robot, config, "ee");
            transformResultOur = obj.directKinematics(q_values(1),q_values(2),q_values(3));
            
            disp("Toolkit Direct Kinematics: ");
            disp(transformResultToolKit);
            
            disp("Calculated Direct Kinematics: ");
            disp(transformResultOur);
            
            disp("Inverse Kinematics: ");
            disp("      [q1,    q2,     q3]: ");
            [q1_ik, q2_ik, q3_ik] = obj.inverseKinematics(q_values(1),q_values(2),q_values(3));
            disp(double(subs([q1_ik, q2_ik, q3_ik], parameters_old,parameters_new)));
            disp(q1_ik); disp(q2_ik); disp(q3_ik);
                        
            disp("Geometrical Jacobian Toolkit: ");
            disp(geometricJacobian(robot, config,"ee"));
            
            disp("Geometrical Jacobian Calcuated: ");
            disp(obj.GeometricalJacobianMatrix)
            disp(obj.solveGeometricalJacobianMatrix(q_values(1),q_values(2),q_values(3)));
            
            disp("Analytical Jacobian Calculated (3x3): ");
            disp(obj.solveAnalyticalJacobianMatrix(q_values(1),q_values(2),q_values(3)));
           
            
            disp("Analytical Jacobian Calculated: (6x3): ");
            disp(obj.AnalyticalJacobianMatrixComplete)
            disp(obj.solveAnalyticalJacobianComplete(q_values(1),q_values(2),q_values(3)));
            
            disp("<<-- Testing Homogenous Transformations -->>")
            
            H1_1_ours = obj.getHomogeneousMatrix(obj.DH_table, 1,1);
            H1_1_tool = getTransform(robot, config, "base_link");
            disp("Toolbox H_1_1 results: ");
            disp(H1_1_tool);
            disp("Calculated H_1_1 results: ");
            disp(double(subs(H1_1_ours, parameters_old,parameters_new)));
            
            H1_2_ours = obj.getHomogeneousMatrix(obj.DH_table, 2,3);
            H1_2_handMade = obj.dhMatrixGenerator(obj.DH_table, 2);
            disp("Handmade H_2_3 results: ");
            disp(double(subs(H1_2_handMade, parameters_old,parameters_new)));
            disp("Calculated H_2_3 results: ");
            disp(double(subs(H1_2_ours, parameters_old,parameters_new)));
            
            H1_3_ours = obj.getHomogeneousMatrix(obj.DH_table, 1,3);
            H1_3_handMade = obj.dhMatrixGenerator(obj.DH_table, 1) * obj.dhMatrixGenerator(obj.DH_table, 2);
            disp("Handmade H_1_3 results: ");
            disp(double(subs(H1_3_handMade, parameters_old,parameters_new)));
            disp("Calculated H_1_3 results: ");
            disp(double(subs(H1_3_ours, parameters_old,parameters_new)));
            
            H4_5_ours = obj.getHomogeneousMatrix(obj.DH_table, 4,5);
            H4_5_handMade = obj.dhMatrixGenerator(obj.DH_table, 4) * [
                [0,0,1,0];
                [0,1,0,0];
                [-1,0,0,0];
                [0,0,0,1];
            ];
            disp("Handmade H_4_5 results: ");
            disp(double(subs(H4_5_handMade, parameters_old,parameters_new)));
            disp("Calculated H_4_5 results: ");
            disp(double(subs(H4_5_ours, parameters_old,parameters_new)));
            
            H_toolbox_ee = getTransform(robot, config, "ee");
            H_ours_ee = obj.getHomogeneousMatrix(obj.DH_table,1,5);
            
            disp("Toolbox H_1_5 (ee) results: ");
            disp(double(subs(H_toolbox_ee, parameters_old,parameters_new)));
            disp("Calculated H_1_5(ee) results: ");
            disp(H_ours_ee);
            disp(double(subs(H_ours_ee, parameters_old,parameters_new)));
            
            disp("<<-- Testing Partials Jacobians -->>")
            disp("J_P_2_3 Handmade: ");
            H = obj.getHomogeneousMatrix(obj.DH_table, 2,3);
            P_C_2_3 = H(1:3,1:3) * obj.P_COM(1,:)' + H(1:3,4);
            
            H_pre = obj.getHomogeneousMatrix(obj.DH_table,2,2); %H_j-1
            P_J_m_1 = H_pre(1:3,4); %P_j-1
            Z2 = H_pre(1:3,3);
            partial = sym(zeros(6,3));
            partial(4:6,1) = cross(Z2, P_C_2_3 - P_J_m_1);
            partial(1:3,1) = Z2;
            disp(simplify(partial));
            
            disp("J_P_2_3 Calculated");
            partial_ours = obj.partialJacobian(obj.DH_table,obj.P_COM,obj.Joints, 2);
            disp(simplify(partial_ours));
            
            disp("J_P_2_4 Handmade: ");
            partial = sym(zeros(6,3));
            
            H = obj.getHomogeneousMatrix(obj.DH_table, 2,4);
            P_C_2_4 = H(1:3,1:3) * obj.P_COM(2,:)' + H(1:3,4);
            
            H_2_2 = obj.getHomogeneousMatrix(obj.DH_table,2,2); %H_j-1
            P_2_2 = H_2_2(1:3,4); %P_j-1
            Z2 = H_pre(1:3,3);
                     
            partial(4:6,1) = cross(Z2, P_C_2_4 - P_2_2);
            partial(1:3,1) = Z2;
            
            H_2_3 = obj.getHomogeneousMatrix(obj.DH_table,2,3); %H_j-1
            partial(4:6,2) = H_2_3(1:3,3);
            
            disp(simplify(partial));
            
            disp("J_P_2_4 Calculated");
            partial_ours = obj.partialJacobian(obj.DH_table,obj.P_COM,obj.Joints, 3);
            disp(simplify(partial_ours));
                        
            disp("J_P_2_5 Handmade: ");
            partial = sym(zeros(6,3));
            
            H = obj.getHomogeneousMatrix(obj.DH_table, 2,5);
            P_C_2_5 = H(1:3,1:3) * obj.P_COM(3,:)' + H(1:3,4);
            
            H_2_2 = obj.getHomogeneousMatrix(obj.DH_table,2,2); %H_j-1
            P_2_2 = H_2_2(1:3,4); %P_j-1
            Z2 = H_pre(1:3,3);
                        
            partial(4:6,1) = cross(Z2, P_C_2_5 - P_2_2);
            partial(1:3,1) = Z2;
            
            H_2_3 = obj.getHomogeneousMatrix(obj.DH_table,2,3); %H_j-1
            partial(4:6,2) = H_2_3(1:3,3);
            
            H_2_4 = obj.getHomogeneousMatrix(obj.DH_table,2,4); %H_j-1
            partial(4:6,3) = H_2_4(1:3,3);
                       
            disp(simplify(partial));
            
            disp("J_P_2_5 Calculated");
            partial_ours = obj.partialJacobian(obj.DH_table,obj.P_COM,obj.Joints, 4);
            disp(simplify(partial_ours));
            
            disp("<<-- Potential Energy -->>");
            PE = obj.potentialEnergy(obj.masses, obj.DH_table,obj.P_COM ,obj.g);
            disp(PE);
            disp(double(subs(PE, parameters_old, parameters_new)));
            
            disp(" <<-- Gravity Matrix -->>");
            disp("Gravity Differential Matrix");
            gravityDifferential = obj.gravityMatrixDrifferential(PE, obj.q);
            disp(gravityDifferential);
            disp(double(subs(gravityDifferential, parameters_old, parameters_new)));
            
            disp("Gravity Matrix");
            gravityPotential = obj.gravityMatrix(obj.DH_table, obj.P_COM, obj.Joints, obj.masses, obj.g);
            disp(gravityPotential)
            disp(double(subs(gravityPotential,  parameters_old, parameters_new)));
            
            disp("<<-- Inertia Matrix -->>");
            Bq = obj.inertiaMatrix(obj.DH_table, obj.masses, obj.P_COM, obj.Joints, obj.lengths);
            Bq_solve = double(subs(Bq, parameters_old, parameters_new));
            disp(Bq)
            disp(Bq_solve);
                   
            if(Bq_solve == Bq_solve')
               disp("Bq matrix is skew-symetric!");
            else
               disp("Bq matrix is not skew-symetric!");
            end
            
            eig_vals = eig(Bq_solve);
            disp("Eigen Values:");
            disp(eig_vals);
            
            is_pd = [sum([sign(eig_vals)==1])==length(Bq_solve)];
            if is_pd ==1 
                disp("Bq is Positive Definite")
            else 
                disp("Bq isn't Positive Definite")
            end
            
            disp("Corolios");
            C = obj.coriolisMatrix(Bq, obj.q, obj.q_dot);
            disp(C)
            % disp(double(subs(C,  parameters_old, parameters_new)));

            disp(" <<-- Kinetic Energy -->>");
            KE = obj.kineticEnergy(Bq, obj.q_dot);
            disp(KE);
            disp(double(subs(KE, parameters_old, parameters_new)));
            
            disp("<<-- Langrangian Equation -->>");
            t_lang = obj.equationOfMotionLangrangian();
            disp(t_lang);
            disp(double(subs(t_lang,parameters_old, parameters_new)));
            
            disp("<<-- Newton Euler -->>");
            t_NE = obj.equationOfMotionNewtonEuler();
            disp(t_NE);
            disp(double(subs(t_NE,parameters_old, parameters_new)));

        end

        function toBeDeleted(obj)
            s = simplify(obj.getHomogeneousMatrix(obj.DH_table,1,5))
            disp("Latex");
            latex(s)
        end
        function EoM = equationOfMotionLangrangian(obj)
            Bq = obj.inertiaMatrix(obj.DH_table, obj.masses, obj.P_COM, obj.Joints, obj.lengths);
            C = obj.coriolisMatrix(Bq, obj.q, obj.q_dot);
            G = obj.gravityMatrix(obj.DH_table, obj.P_COM, obj.Joints, obj.masses, obj.g);

            EoM = Bq*obj.q_dot_dot + C*obj.q_dot + diag(obj.Fv)*obj.q_dot + diag(obj.Fs)*sign(obj.q_dot) + G;
        end

        function EoM = equationOfMotionNewtonEuler(obj)
            [w,w_dot, p_dot_dot_c] = obj.forwardEquation(obj.g, obj.DH_table, obj.q_dot, obj.q_dot_dot, obj.P_COM, obj.Joints);
            EoM = obj.backwordEquation(w,w_dot,p_dot_dot_c, obj.f_e, obj.mu_e, obj.DH_table, obj.P_COM, obj.masses,obj.Joints, obj.lengths, obj.Fv, obj.Fs, obj.q_dot)';
        end

        function res = solveAnalyticalJacobianComplete(obj, t_1, d_1, d_2)
            analyticalJacobian = obj.AnalyticalJacobianMatrixComplete;
            x = obj.TransoformationMatrix;

            rotMatrix = double(subs(x, [obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3),obj.q(1), obj.q(3),obj.q(2)], [0.4, 0.4, 0.46, 0.15, t_1,d_2, d_1]));
            phi2 = rotm2eul(rotMatrix(1:3,1:3),"ZYZ");
            
            % res = analyticalJacobian;
            res = double(subs(analyticalJacobian, [obj.phi(1),obj.phi(2),obj.phi(3),obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3),obj.q(1), obj.q(3),obj.q(2)], [phi2(1),phi2(2),phi2(3),0.15, 0.4, 0.3, 0.16, t_1,d_2, d_1]));

        end

        function testOperationalSpace(obj, q_values)

           density = 2700; %Almunium density
           mass_values = [pi * obj.lengths(1,1)^2*obj.lengths(1,3)*density, obj.lengths(2,1)*obj.lengths(2,2)*obj.lengths(2,3)*density, obj.lengths(3,1)*obj.lengths(3,2)*obj.lengths(3,3)*density];
           lenghts_values = [
               [0.02, 2*pi*0.02, 0.4]; 
               [0.03, 0.03, 0.46]; 
               [0.01, 0.01, 0.15]
               ];
           mass_values = subs(mass_values, [obj.lengths(1,1),obj.lengths(1,2), obj.lengths(1,3), obj.lengths(2,1),obj.lengths(2,2), obj.lengths(2,3),obj.lengths(3,1),obj.lengths(3,2), obj.lengths(3,3)], [lenghts_values(1,1), lenghts_values(1,2), lenghts_values(1,3),lenghts_values(2,1), lenghts_values(2,2), lenghts_values(2,3),lenghts_values(3,1), lenghts_values(3,2), lenghts_values(3,3)]);
           q_dot_values = [0;0;0];
           q_dot_dot_values = [0;0;0];
           g_value = 9.806;
           f_e_values = [0;0;0];
           mu_e_values = [0;0;0];
            
           x = obj.getHomogeneousMatrix(obj.DH_table, 2, 5);

           rotMatrix = double(subs(x, [obj.r1,obj.lengths(1,3),obj.lengths(2,3),obj.lengths(3,3),obj.q(1), obj.q(3),obj.q(2)], [0.4, 0.4, 0.46, 0.15, q_values(1),q_values(2), q_values(3)]));
           phi_vals = rotm2eul(rotMatrix(1:3,1:3),"ZYZ");

           x_dot_dot_valyes = [0;0;0;0;0;0];

           parameters_old = [obj.lengths(1,1),obj.lengths(1,2), obj.lengths(1,3), obj.lengths(2,1),obj.lengths(2,2), obj.lengths(2,3),obj.lengths(3,1),obj.lengths(3,2), obj.lengths(3,3), obj.masses, obj.g, obj.q', obj.q_dot', obj.q_dot_dot', obj.phi, obj.x_dot_dot'];
           parameters_new  = [lenghts_values(1,1), lenghts_values(1,2), lenghts_values(1,3),lenghts_values(2,1), lenghts_values(2,2), lenghts_values(2,3),lenghts_values(3,1), lenghts_values(3,2), lenghts_values(3,3), mass_values, g_value, q_values, q_dot_values', q_dot_dot_values',phi_vals, x_dot_dot_valyes'];

           Bq = obj.inertiaMatrix(obj.DH_table, obj.masses, obj.P_COM, obj.Joints, obj.lengths);
           C = obj.coriolisMatrix(Bq, obj.q, obj.q_dot);
           G = obj.gravityMatrix(obj.DH_table, obj.P_COM, obj.Joints, obj.masses, obj.g);

           GeoJacob = obj.generateGeometricalJacobianFromFirstJoint(obj.DH_table, obj.Joints, 4);
           [Ta, AnalJacob] = obj.generateAnalyticalJacobianComplete(GeoJacob);
           [Ba, Ca, Ga] = obj.operationalSpaceDynamicModel(AnalJacob, Bq,C,G, obj.q, obj.q_dot);

           Ba_solved = double(subs(Ba, parameters_old, parameters_new));
           Ca_solved = double(subs(Ca, parameters_old, parameters_new));
           Ga_solved = double(subs(Ga, parameters_old, parameters_new));

            ya = Ba*obj.x_dot_dot + Ca + Ga;

            ya_solved = double(subs(ya,parameters_old,parameters_new));

            
            J_fixed = sym(zeros(6,3));
            
            J_fixed(1:3,:) = GeoJacob(4:6,:);
            J_fixed(4:6,:) = GeoJacob(1:3,:);


            Tau = obj.equationOfMotionLangrangian;

            ya_expected = (Ta'*pinv(J_fixed')*Tau);

            ya_ex_solved = double(subs(ya_expected, parameters_old,parameters_new));

            disp("<<-- OperationalSpace -->>");
            disp("Ba:");
            disp(Ba_solved);
            disp(Ba);
            disp("Ca");
            disp(Ca_solved);
            disp(Ca);
            disp("Ga:");
            disp(Ga_solved);
            disp(Ga);
            disp("ya_expected:");
            disp(double(subs(ya_expected,parameters_old, parameters_new)));
            disp("ya_solved:")
            disp(double(subs(ya_solved,parameters_old, parameters_new)));
        end
    end
    methods (Access = private)

        function matrix = dhMatrixGenerator(~, dh_table, i)
            theta = dh_table(i,1);
            d = dh_table(i,4);
            a = dh_table(i,3);
            alpha = dh_table(i,2);
            
            matrix = [
                [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta)];
                [sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta)];
                [0,          sin(alpha),             cos(alpha),             d];
                [0,          0,                      0,                      1];
            ];

          end

        %Matrix from 1 to index, 1 is base frame
        function matrix = getTransformationMatrix(obj, dh_table, index)
            matrix = eye(4);
            for i = 1:index
                matrix = matrix * obj.dhMatrixGenerator(dh_table, i);
            end
        end

        function jacobianMatrix = generateJacobianMatrix(obj, dh_table,joints, index)
            jacobianMatrix = sym(zeros(6,size(joints,2)));
            pn = obj.getTransformationMatrix(dh_table,index+1);
            pn =  pn(1:3,4); %extracting only position part
            matrix = obj.dhMatrixGenerator(dh_table,1);
            for i = 1:index %1 is base frame
                if joints(i) == "Revolut"
                    jacobianMatrix(1:3,i) = matrix(1:3,3);
                    jacobianMatrix(4:6,i) = cross(matrix(1:3,3),pn - matrix(1:3,4));
                else
                    jacobianMatrix(1:3,i) = [0,0,0];
                    jacobianMatrix(4:6,i) = matrix(1:3,3);
                end
                matrix = matrix * obj.dhMatrixGenerator(dh_table,i+1);
            end

            jacobianMatrix = simplify(jacobianMatrix);
        end
      
        function analyticalJacobianMatrix = generateAnalyticalJacobianMatrix(~,pe,q)
            analyticalJacobianMatrix = sym(zeros(3,size(q,1)));
            
            for j = 1:3
               for i = 1:size(q,1)
                   analyticalJacobianMatrix(j,i) = diff(pe(j),q(i));
               end
            end

            analyticalJacobianMatrix = simplify(analyticalJacobianMatrix);
        end

        %1 is the base frame
        function H = getHomogeneousMatrix(obj, dh_table, startIndex, endIndex)
            H = sym(eye(4,4));

            if(startIndex == endIndex) %if no transformation is needed, use empty homogeneous transformation
                return;
            end

            if(endIndex > 5)
                return;
            end

            for i = startIndex: endIndex-1
                H = H * obj.dhMatrixGenerator(dh_table, i);
            end

            H = simplify(H);
        end

        %index 1 is base, but in this the base is transfered to link 1 so
        %it is 2
        function partial = partialJacobian(obj, dh_table, P_COM, joints, index)
            partial = sym(zeros(6, size(joints,2)));
            H = obj.getHomogeneousMatrix(dh_table, 2, index+1);
            P_li = H(1:3,1:3)*P_COM(index-1,:)'+H(1:3,4);

            for i=2:index
                H_pre = obj.getHomogeneousMatrix(dh_table,2, i);

                if(joints(i-1) == "Revolut")
                    partial(4:6,i-1) = cross(H_pre(1:3,3), P_li - H_pre(1:3,4));
                    partial(1:3,i-1) = H_pre(1:3,3);
                else
                    partial(4:6,i-1) = H_pre(1:3,3);
                end

            end
            partial = simplify(partial);
        end

        %index 1 is base, but in this the base is transfered to link 1 so
        %it is 2
        function partial = generateGeometricalJacobianFromFirstJoint(obj, dh_table, joints, index)
            partial = sym(zeros(6, size(joints,2)));
            H = obj.getHomogeneousMatrix(dh_table, 2, index+1);
            P_li = H(1:3,4);

            for i=2:index
                H_pre = obj.getHomogeneousMatrix(dh_table,2, i);

                if(joints(i-1) == "Revolut")
                    partial(4:6,i-1) = cross(H_pre(1:3,3), P_li - H_pre(1:3,4));
                    partial(1:3,i-1) = H_pre(1:3,3);
                else
                    partial(4:6,i-1) = H_pre(1:3,3);
                end

            end
            partial = simplify(partial);
        end
    
        function PE = potentialEnergy(obj, masses, dh_table, CoM, g)
            PE = 0;
            gvec = [0;-g;0];
            for i = 2:4
                H = obj.getHomogeneousMatrix(dh_table, 2,i+1);
                Pl = H(1:3,1:3)*CoM(i-1,:)'+H(1:3,4);
                PE = PE + (-masses(i-1)*gvec'*Pl);
            end
            
            PE = simplify(PE);
        end

        function G = gravityMatrixDrifferential(~, U, q)
            G = [diff(U,q(1)); diff(U,q(2)); diff(U,q(3))];
            G = simplify(G);
        end

        function G = gravityMatrix(obj, dh_table, P_COM, joints, masses, g)
            G = sym(zeros(3,1));
            gvec = [0;-g;0];
            for i = 2:4
                for j = 2:4
                    Jp = obj.partialJacobian(dh_table,P_COM, joints, j);
                    G(i-1,1) = G(i-1,1) + (-masses(j-1)*gvec'*Jp(4:6, i-1));
                end
            end
        end
    
        function Bq = inertiaMatrix(obj, dh_table, masses, P_COM, joints, lengths)
            Bq = 0;
            for i=2:4
                H = obj.getHomogeneousMatrix(dh_table, 2, i+1);
                J = obj.partialJacobian(dh_table, P_COM, joints, i);
                if(joints(i-1) == "Revolut")
                    Ic = obj.inertiaCylender(lengths(i-1,1),lengths(i-1,2),lengths(i-1,3),masses(i-1));
                else
                    Ic = obj.inertiaCube(lengths(i-1,1),lengths(i-1,2),lengths(i-1,3),masses(i-1));
                end
                I = obj.translateInertia(Ic, masses(i-1),P_COM(i-1,:)');
                Bq = Bq + (masses(i-1)*J(4:6,:)'*J(4:6,:)+ J(1:3,:)'*H(1:3,1:3)*I*H(1:3,1:3)'*J(1:3,:));
            end
            
            Bq = simplify(Bq);
        end

        function Ic = inertiaCylender(obj, a,b,c,m)
            Ic = [
                [1/2*m*(a^2+b^2),0,0];
                [0,1/2*m*(3*(a^2+b^2)^2 + c^2),0];
                [0,0,1/2*m*(3*(a^2+b^2)^2 + c^2)]
            ];
        end

        function Ic = inertiaCube(obj, a,b,c,m)
            Ic = [
                [1/12*m*(b^2+a^2),0,0];
                [0,1/12*m*(a^2+c^2),0];
                [0,0,1/12*m*(c^2+b^2)]
            ];
        end

        function I = translateInertia(obj, Ic, m,r)
            I = Ic + m*(r'*r*eye(3,3)-r*r');
        end
        
        function KE = kineticEnergy(obj, Bq, q_dot)
            KE = 1/2*q_dot'*Bq*q_dot;
        end
    
        function C = coriolisMatrix(obj, Bq, q, q_dot)
            C = sym(zeros(3,3));
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        C(i,j) = C(i,j) + (1/2 *(diff(Bq(i,j),q(k))+diff(Bq(i,k),q(j))-diff(Bq(j,k),q(i)))*q_dot(k));
                    end
                end
            end
            C = simplify(C);
        end
    
        function [w_app, w_dot_app, p_dot_dot_c_app] = forwardEquation(obj, g, dh_table, q_dot, q_dot_dot,P_COM, joints)
            w_pre = [0;0;0];
            w_dot_pre = [0;0;0];
            p_dot_dot_pre = [0;g;0];
            
            Z0 = [0;0;1];

            p_dot_dot_c_app = sym(zeros(3,1,3));
            w_app = sym(zeros(3,1,3));
            w_dot_app = sym(zeros(3,1,3));

            for i = 2:4
                H = obj.getHomogeneousMatrix(dh_table, i,i+1);
                R = H(1:3,1:3);
                w = R' * w_pre;
                w_dot = R' * w_dot_pre;
                if(joints(i-1) == "Revolut")
                    w = w + R'*q_dot(i-1)*Z0;
                    w_dot = w_dot + R' * (q_dot_dot(i-1)*Z0 + cross(q_dot(i-1)*w_pre, Z0));
                end
                rlink = R' * H(1:3,4);
                p_dot_dot = R' * p_dot_dot_pre + cross(w_dot, rlink) + cross(w,cross(w,rlink));
                if(joints(i-1) == "Prismatic")
                    p_dot_dot = p_dot_dot + R' * q_dot_dot(i-1)*Z0 + cross(2*q_dot(i-1)*w,R'*Z0);
                end
                p_dot_dot_c = p_dot_dot + cross(w_dot,P_COM(i-1,1:3)') + cross(w,cross(w,P_COM(i-1,1:3)'));

                %appending the results
                p_dot_dot_c_app(:,:,i-1) = p_dot_dot_c;
                w_dot_app(:,:,i-1) = w_dot;
                w_app(:,:,i-1) = w;

                %replacing old with new
                w_pre = w;
                w_dot_pre = w_dot;
                p_dot_dot_pre = p_dot_dot;
            end
        end
        
        function tourqe = backwordEquation(obj, w,w_dot,p_dot_dot_c, f_e, mu_e, dh_table, P_COM, masses, joints, lengths, Fv, Fs, q_dot)
            f_next = f_e;
            mu_next = mu_e;
            Z0 = [0;0;1];

            for i = 5:-1:3
                H = obj.getHomogeneousMatrix(dh_table, i,i+1);
                R = H(1:3,1:3);

                H_m_1 = obj.getHomogeneousMatrix(dh_table, i-1,i);
                R_m_1 = H_m_1(1:3,1:3);

                rlink = R_m_1' * H_m_1(1:3,4);
                f = R*f_next + masses(i-2) * p_dot_dot_c(:,:,i-2);

                if(joints(i-2) == "Revolut")
                    Ic = obj.inertiaCylender(lengths(i-2,1),lengths(i-2,2),lengths(i-2,3),masses(i-2));
                else
                    Ic = obj.inertiaCube(lengths(i-2,1),lengths(i-2,2),lengths(i-2,3),masses(i-2));
                end
                I = obj.translateInertia(Ic, masses(i-2),P_COM(i-2,:)');

                mu = -cross(f,rlink + P_COM(i-2,:)') + R*mu_next + cross(R*f_next,P_COM(i-2,:)') + I*w_dot(:,:,i-2) + cross(w(:,:,i-2),I*w(:,:,i-2));

                if(joints(i-2) == "Prismatic")
                    tourqe(i-2) = f' * R_m_1'*Z0;
                else
                    tourqe(i-2) = mu' * R_m_1'*Z0;
                end

                tourqe(i-2) = tourqe(i-2) + Fv(i-2)*q_dot(i-2) + Fs(i-2)*sign(q_dot(i-2));

                f_next = f;
                mu_next = mu;
            end
        end
    
        function [Ta, analyticalJacobian] = generateAnalyticalJacobianComplete(obj, J)
            
            T = [[0, -sin(obj.phi(1)), cos(obj.phi(1))*sin(obj.phi(2))];
                 [0, cos(obj.phi(1)), sin(obj.phi(1))*sin(obj.phi(2))];
                 [1,0,cos(obj.phi(2))]];

            Ta = simplify([eye(3) zeros(3);
                  zeros(3) T]);
            Jac = sym(zeros(6,3));
            Jac(1:3,1:3) = J(4:6, 1:3);
            Jac(4:6,1:3) = J(1:3,1:3);

            res  = simplify(Ta \ Jac);

            analyticalJacobian = sym(zeros(6,3));

            analyticalJacobian(1:3,:) = res(4:6,:);
            analyticalJacobian(4:6,:) = res(1:3,:);

            Ta = simplify(Ta);
        end
        
        function [Ba, Ca, Ga] = operationalSpaceDynamicModel(obj, Ja, B, C, G, q, q_dot)
            analyticalJacobian = sym(zeros(6,3));
            analyticalJacobian(1:3,:) = Ja(4:6,:);
            analyticalJacobian(4:6,:) = Ja(1:3,:);

            Ja_inv = pinv(analyticalJacobian);

            time = sym("t");
            q_time = [symfun("q1(t)", time); symfun("q2(t)", time); symfun("q3(t)", time)];
            q_time_dot = diff(q_time);

            J_dot = subs(analyticalJacobian, q, q_time);
            J_dot = diff(J_dot);

            J_dot = subs(J_dot, [q_time_dot; q_time], [q_dot; q]);

            Ba = Ja_inv' * B * Ja_inv;
            Ca = Ja_inv' * C * q_dot - Ba *J_dot*q_dot;
            Ga = Ja_inv' * G;
        end
    end
end
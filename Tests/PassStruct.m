function f = PassStruct(a, i, phi, mesh, matrices, params, model)

model.B_mu.a_w(i) = a;
f = Valve_GetJ(phi, mesh, matrices, params, model);

end


%TODO: if nargout >= 2
%     
%     compute derivative wrt a
%     
% end
% 
% [a,b] = qwe();
% a = qwe();


"""
Rotation Matricies
    input angles defined in radians
"""
function rotz(γ)
    rotmat = [cos(γ) -sin(γ) 0; sin(γ) cos(γ) 0; 0 0 1];
    return rotmat
end
function rotx(α)
    rotmat = [1 0 0;0 cos(α) -sin(α); 0 sin(α) cos(α)];
    return rotmat
end

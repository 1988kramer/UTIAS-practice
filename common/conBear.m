function newBear = conBear(oldBear)
    while oldBear < -pi
        oldBear = oldBear + 2*pi;
    end
    while oldBear > pi
        oldBear = oldBear - 2*pi;
    end
    newBear = oldBear;
end
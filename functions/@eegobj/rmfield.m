function this = rmfield(this, vals);

    eval( [ 'this.EEG.' vals '=[];' ] );

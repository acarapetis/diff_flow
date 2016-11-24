var grid_size = 10;
var draw_scale = 50;
var canvas = $('canvas')[0];
var ctx = canvas.getContext('2d');
var real_scale = draw_scale;
var resize = function() {
    // Make sure 1 canvas pixel = 1 screen pixel
    var dpr = window.devicePixelRatio || 1;
    var bsr = ctx.webkitBackingStorePixelRatio 
          ||  ctx.mozBackingStorePixelRatio
          ||  ctx.msBackingStorePixelRatio
          ||  ctx.oBackingStorePixelRatio
          ||  ctx.backingStorePixelRatio || 1;
    var PIXEL_RATIO = dpr/bsr;

    canvas.width    = $('canvas').width() * PIXEL_RATIO;
    canvas.height   = $('canvas').height() * PIXEL_RATIO;
    ctx.translate(0.5,0.5);
    real_scale = PIXEL_RATIO * draw_scale;
};

// Fix mod to target [0,n)
var mod = function(m,n) {
    return ((m%n)+n)%n;
};

var Torus = function (w,h,initf) {
    this.width = w;
    this.height = h;

    this.data = [];

    if (initf === undefined) {
        initf = function(x,y) { return null; };
    }

    for (var i = 0; i < w; i++) {
        var col = [];
        for (var j = 0; j < h; j++) {
            col.push(initf(i,j));
        }
        this.data.push(col);
    }
};

Torus.prototype.el = function (x,y) {
    x = mod(x,this.width);
    y = mod(y,this.height);
    return this.data[x][y];
};

Torus.prototype.each_el = function(f) {
    for (var i = 0; i < this.width; i++) {
        for (var j = 0; j < this.width; j++) {
            f(this.data[i][j],i,j);
        }
    }
}

Torus.prototype.mutate_el = function(f) {
    for (var i = 0; i < this.width; i++) {
        for (var j = 0; j < this.width; j++) {
            this.data[i][j] = f(this.data[i][j],i,j);
        }
    }
}

Torus.prototype.neighbours = function (x,y) {
    return [
        this.el(x,y+1),
        this.el(x,y-1),
        this.el(x+1,y),
        this.el(x-1,y)
    ];
};

Torus.prototype.SEneighbours = function(x,y) {
    return [
        this.el(x,y+1),
        this.el(x+1,y),
    ];
};

Torus.prototype.grad = function (x,y,f) {
    // Central difference approximation to grad(f) @ (x,y)
    return [
        (f(this.el(x+1,y)) - f(this.el(x-1,y)))/2,
        (f(this.el(x,y+1)) - f(this.el(x,y+1)))/2
    ];
};

var xc = function(e) { return e.x; };
var yc = function(e) { return e.y; };

var vadd = function(u,v) {
    return {x: u.x + v.x, y: u.y + v.y};
};

var vscale = function(u,k) {
    return {x: k*u.x , y: k*u.y};
};

var s1diff = function(a,b,l) {
    // The "displacement" a-b on a circle of length l
    return mod((mod(a,l)-mod(b,l))+l/2,l)-l/2;
};

Torus.prototype.minus = function(u,v) {
    return {x: s1diff(u.x,v.x,this.width), y: s1diff(u.y,v.y,this.height)};
};

Torus.prototype.d2 = function(u,v) {
    var d = this.minus(u,v);
    return d.x*d.x + d.y*d.y;
}


Torus.prototype.jac = function(x,y) {
    // Central difference approximation to Jacobian matrix
    var e0 = this.el(x,y);
    var tt = this;
    return [ 
        this.grad(x,y,function(e) { return s1diff(e.x,e0.x,tt.width); }), 
        this.grad(x,y,function(e) { return s1diff(e.y,e0.y,tt.height); })
    ];
};

Torus.prototype.lap = function(x,y) {
    // 2nd difference approximation to harmonic map Laplacian
    var e0 = this.el(x,y);
    return this.neighbours(x,y)
               .map(function(e){ return this.minus(e,e0); }.bind(this))
               .reduce(function(a,b){ return vadd(a,b); },{x:0,y:0});
};

Torus.prototype.draw_edge_on = function(c,u,e) {
    var d = tt.minus(u,e);
    var z = vadd(e,d);
    var do_edge = function(x1,y1,x2,y2) {
        c.beginPath();
        c.moveTo(x1 * real_scale, y1 * real_scale);
        c.lineTo(x2 * real_scale, y2 * real_scale);
        c.stroke();
    };
    
    // Draw fundamental representative:
    do_edge(e.x,e.y,z.x,z.y);

    var w = this.width;
    var h = this.height;

    // Now draw the necessary translates.
    // First the Von Neumann neighbourhood:
    if (z.x > w)
        do_edge(e.x-w, e.y, z.x-w, z.y);

    if (z.x < 0)
        do_edge(e.x+w, e.y, z.x+w, z.y);

    if (z.y > h)
        do_edge(e.x, e.y-h, z.x, z.y-h);

    if (z.y < 0)
        do_edge(e.x, e.y+h, z.x, z.y+h);

    // Finally the diagonals:
    if (z.x > w && z.y > h)
        do_edge(e.x-w, e.y-h, z.x-w, z.y-h);

    if (z.x > w && z.y < 0)
        do_edge(e.x-w, e.y+h, z.x-w, z.y+h);

    if (z.x < 0 && z.y < 0)
        do_edge(e.x+w, e.y+h, z.x+w, z.y+h);

    if (z.x < 0 && z.y > h)
        do_edge(e.x+w, e.y-h, z.x+w, z.y-h);

};

Torus.prototype.draw_on = function(c) {
    tt = this;
    this.each_el(function(e,x,y) {
        tt.SEneighbours(x,y).forEach(function(e2) {
            c.strokeStyle='#000000';
            c.lineWidth=1;
            tt.draw_edge_on(c,e2,e);
        });
    });
        /*
    this.each_el(function(e,x,y) {
        var v = tt.lap(x,y);
        
        c.beginPath();
        c.arc(e.x * real_scale,e.y * real_scale, 0.1 * real_scale, 0, 2*Math.PI);
        c.closePath();
        c.fill();

        c.strokeStyle='#FF0000';
        c.lineWidth=4;
        c.beginPath();
        c.moveTo(e.x * real_scale, e.y * real_scale);
        c.lineTo((e.x + 10*v.x) * real_scale, (e.y+10*v.y)*real_scale);
        c.stroke();
    });*/
};

Torus.prototype.crop = function(e) {
    return { x: mod(e.x, this.width), y: mod(e.y, this.height) };
}

var tor_add = function(t,s) {
    return new Torus(t.width,t.height,function(x,y) { 
        return vadd(t.el(x,y), s.el(x,y));
    });
}

var hmhf_step = function(t, dt) {
    return new Torus(t.width, t.height, function(x,y) {
        return vadd(t.el(x,y), vscale(t.lap(x,y),dt));
    });
}

var myflow_step = function(t, dt) {
    return new Torus(t.width, t.height, function(x,y) {
        var j = t.jac(x,y);
        var r = j[0][0] * j[0][0] + j[0][1] * j[0][1] + j[1][0] * j[1][0] + j[1][1] * j[1][1]
                + 2*(j[0][0]*j[1][1] - j[0][1]*j[1][0]);
        return t.crop(vadd(t.el(x,y), vscale(t.lap(x,y),dt/r)));
    });
}

/*
ctx.tdraw_circle = function(x,y,r) {
    x = mod(x,
    c.beginPath();
    c.arc(x * draw_scale,y * draw_scale, r * draw_scale, 0, 2*Math.PI);
    c.closePath();
    c.fill();
}
*/

var dragging = false;
var drag_origin;
var mouse_pos;

var bump = function(x2) {
    return x2 > 1 ? 0 : Math.exp(-1/(1-x2));
}

var mousemove = function(evt) {
    var rect = canvas.getBoundingClientRect();
    var old_pos = mouse_pos;
    mouse_pos = {
        x: (evt.clientX-rect.left)/(rect.right-rect.left)*canvas.width/real_scale,
        y: (evt.clientY-rect.top)/(rect.bottom-rect.top)*canvas.height/real_scale
    };
    if (dragging) {
        var impulse = t.minus(mouse_pos, old_pos);
        var drag_coeff = $('#dcoeff')[0].value;
        t.mutate_el(function(e1,i,j) {
            var scal = bump(t.d2(e1,old_pos)/(0.01*t.width*t.width*drag_coeff*drag_coeff))
            return t.crop(vadd(e1, vscale(impulse,scal)));
        });
    }
};
$(canvas).on('mousemove',mousemove);
$(canvas).on('touchmove', function(e) { 
    mousemove(e.originalEvent.changedTouches[0]); 
    return false;
});

var mousedown = function(e) {
    if (('button' in e) && e.button > 0) return;
    mousemove(e);
    dragging = true;
    drag_origin = mouse_pos;
    return false;
}
$(canvas).on('mousedown',mousedown);
$(canvas).on('touchstart', function(e) {
    mousedown(e.originalEvent.changedTouches[0]); 
    return false;
});

$(canvas).on('mouseup touchend', function(e) {
    if (!dragging) return;
    if (('button' in e) && e.button > 0) return;
    dragging = false;
    /*
    var drag_terminus = mouse_pos;
    var impulse = t.minus(drag_terminus, drag_origin);
    var drag_coeff = $('#dcoeff')[0].value;
    t.mutate_el(function(e1,i,j) {
        var scal = Math.exp(-drag_coeff*t.d2(e1,drag_origin)/(t.width*t.height));
        return t.crop(vadd(e1, vscale(impulse,scal)));
    });
    draw();
    */
});


var t;

var draw = function() {
    ctx.clearRect(0,0,canvas.width,canvas.height);
    t.draw_on(ctx);
};
var tick = function() {
    draw();
    if ($('#toggle_flow')[0].checked) {
        t = myflow_step(t,$('#timestep')[0].value/100);
    }
};
var reset = function(f) {
    console.log(grid_size);
    if (f === undefined) {
        f = function(x,y) { return { x: x, y: y }; };
    }
    t = new Torus(grid_size,grid_size,f);
    draw();
};
$('#reset').on('click', function() { 
    reset(function(x,y){
        return {x: x, y:y};
    });
});

$('#shear').on('click', function() { 
    reset(function(x,y){
        return t.crop({x: x, y:x+y});
    });
});

var resize_reset = function() {
    grid_size = parseInt($('#size')[0].value);
    draw_scale = 300/grid_size;
    resize();
    reset();
};

$('#size').on('input change', resize_reset);

$(window).on('resize orientationchange', function() { resize(); draw(); });
$(window).on('load',function() {
    grid_size = parseInt($('#size')[0].value);
    draw_scale = 300/grid_size;
    $('canvas').width(grid_size * draw_scale).height(grid_size * draw_scale);
    resize_reset();
    t = new Torus(grid_size,grid_size,function(x,y) { return {
        x: x +0.4*Math.cos(2*Math.PI*y/grid_size),
        y: y +0.4*Math.sin(2*Math.PI*x/grid_size)
    }; });
    draw();
    setInterval(tick,15);
});

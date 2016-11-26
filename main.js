// {{{ 2x2 Matrices
var M22 = function(a,b,c,d) {
    this.a = a;
    this.b = b;
    this.c = c;
    this.d = d;
};

M22.prototype.multL = function(o) { // o * this
    return new M22(
            o.a * this.a + o.b * this.c,
            o.a * this.b + o.b * this.d,
            o.c * this.a + o.d * this.c,
            o.c * this.b + o.d * this.d
    );
};

M22.prototype.transform = function(v) { // this * v
    return {
        x: this.a * v.x + this.b * v.y,
        y: this.c * v.x + this.d * v.y
    };
};

M22.prototype.scale = function(s) { // s this
    return new M22(a*s,b*s,c*s,d*s);
};
// }}}
// {{{ (Non-projective) Modular Group
// Generators S,T,-
// S = J, T = shear, - = -Id
// Ss = Tt = -- = Id
// Relations S^2 = (ST)^3 = -
var MS = new M22(0,-1,1,0);
var MSinv = new M22(0,1,-1,0);
var MT = new M22(1,1,0,1);
var MTinv = new M22(1,-1,0,1);
var Mneg = new M22(-1,0,0,-1);
var Mid = new M22(1,0,0,1);

var word_reduce = function(w) {
    // We're not solving the word problem here. Just doing our best.
    var v;
    do {
        v = w;
        w = w.replace(/SS|ss|STSTST|tststs/g,'-');
        w = w.replace(/Ss|sS|Tt|tT/g,'');
        minuses = (w.match(/-/g) || []).length;
        w = (minuses % 2 == 0 ? '' : '-') + w.replace(/-/g,'');
    } while (v != w);
    return w;
};

var Word = function(w) {
    this.w = word_reduce(w) || '';
};

Word.prototype.toString = function() {
    return this.w;
};

var char_to_matrix = function(c) {
    mapping = {
        '-' : Mneg,
        S: MS, s: MSinv,
        T: MT, t: MTinv
    };
    return mapping[c];
}

Word.prototype.toMatrix = function() {
    return this.w.split('').reduceRight(function(r,l) {
        return r.multL(char_to_matrix(l));
    }, Mid);
};
// }}}
var grid_size = 10;
var draw_scale = 50;
var canvas = $('canvas')[0];
var renderer, plane, gfx, texture;
var glw, glh;
var t, idt;
var coloured_grid = false;
var pixel_ratio = function() {
    var dpr = window.devicePixelRatio || 1;
    var ctx;
    try { ctx = canvas.getContext('webgl'); }
    catch(e) { ctx = canvas.getContext('2d'); }

    var bsr = ctx.webkitBackingStorePixelRatio 
          ||  ctx.mozBackingStorePixelRatio
          ||  ctx.msBackingStorePixelRatio
          ||  ctx.oBackingStorePixelRatio
          ||  ctx.backingStorePixelRatio || 1;
    console.log(dpr/bsr);
    return dpr/bsr;
}


// Fix mod to target [0,n)
var mod = function(m,n) {
    return ((m%n)+n)%n;
};

var Torus = function (w,h,initf) {
    this.width = w;
    this.height = h;

    this.data = [];

    if (initf === undefined) {
        initf = function(x,y) { return { x: x, y: y }; };
    }

    for (var i = 0; i < w; i++) {
        var col = [];
        for (var j = 0; j < h; j++) {
            //col.push(this.crop(initf(i,j)));
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

var hsl = function(h, s, l) {
    var r, g, b;

    if(s == 0){
        r = g = b = l; // achromatic
    }else{
        var hue2rgb = function hue2rgb(p, q, t){
            if(t < 0) t += 1;
            if(t > 1) t -= 1;
            if(t < 1/6) return p + (q - p) * 6 * t;
            if(t < 1/2) return q;
            if(t < 2/3) return p + (q - p) * (2/3 - t) * 6;
            return p;
        }

        var q = l < 0.5 ? l * (1 + s) : l + s - l * s;
        var p = 2 * l - q;
        r = hue2rgb(p, q, h + 1/3);
        g = hue2rgb(p, q, h);
        b = hue2rgb(p, q, h - 1/3);
    }

    return 0x10000 * Math.round(r * 255) + 0x100 * Math.round(g * 255)+ Math.round(b * 255);
}

Torus.prototype.domain_colour = function (x,y) {
    return hsl(x/this.width, 1, 0.6 + 0.3*Math.sin(2*Math.PI*y/this.height));
};

Torus.prototype.CM = function() {
    var xb = 0, yb = 0;
    this.each_el(function(e) {
        xb += e.x;
        yb += e.y;
    });
    xb /= this.width * this.height;
    yb /= this.width * this.height;
    return {x: xb, y: yb};
}

Torus.prototype.translateAll = function(dx,dy) {
    this.mutate_el(function(e) {
        return { x: e.x + dx, y: e.y + dy};
    });
}

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
    var d = this.minus(u,e);
    var z = vadd(e,d);
    var do_edge = function(x1,y1,x2,y2) {
        c.moveTo(x1 * real_scale, y1 * real_scale);
        c.lineTo(x2 * real_scale, y2 * real_scale);
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
    c.lineStyle(1, 0, 1);
    this.each_el(function(e,x,y) {
        tt.SEneighbours(x,y).forEach(function(e2) {
            if (coloured_grid) {
                c.lineStyle(2,tt.domain_colour(x,y),1);
            }
            tt.draw_edge_on(c,tt.crop(e2),tt.crop(e));
        });
    });
};

Torus.prototype.crop = function(e) {
    return { x: mod(e.x, this.width), y: mod(e.y, this.height) };
}

Torus.prototype.update_plane = function(p) {
	var t = this;
    for (var x = 0; x <= t.width; x++) {
        for (var y = 0; y <= t.height; y++) {
            var e = t.el(x%t.width,y%t.height);
			var ee = e;
            if (x == t.width || y == t.height) {
                var prev = (x == t.width) ? 
					((y == t.width) ? t.el(t.width - 1, t.height - 1) : t.el(t.width - 1, y%t.height)) : 
					t.el(x%t.height, t.height - 1);
				ee = vadd(prev,t.minus(e,prev));
            }
            p.vertices[2*(t.width+1)*y + 2*x] = (ee.x/t.width)*glw;
            p.vertices[2*(t.width+1)*y + 2*x+1] = (ee.y/t.height)*glh;
        }
    }
};


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
        // r = (u1 + u2)^2
        var r = j[0][0] * j[0][0] + j[0][1] * j[0][1] 
                + j[1][0] * j[1][0] + j[1][1] * j[1][1]
                + 2*(j[0][0]*j[1][1] - j[0][1]*j[1][0]);
        return (vadd(t.el(x,y), vscale(t.lap(x,y),dt/r)));
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

var mouse_pos;

var bump = function(x2) {
    return x2 > 1 ? 0 : Math.exp(-1/(1-x2));
}

var mousemove = function(evt, touch) {
    var rect = canvas.getBoundingClientRect();
    var old_pos = mouse_pos;
    mouse_pos = {
        x: (evt.clientX-rect.left)/(rect.right-rect.left) * t.width,
        y: (evt.clientY-rect.top)/(rect.bottom-rect.top) * t.height
    };
    if (('buttons' in evt && evt.buttons == 1) || (touch && old_pos !== false)) {
        var impulse = t.minus(mouse_pos, old_pos);
        var drag_coeff = $('#dcoeff')[0].value;
        t.mutate_el(function(e1,i,j) {
            var scal = bump(t.d2(e1,old_pos)/(0.01*t.width*t.width*drag_coeff*drag_coeff))
            return vadd(e1, vscale(impulse,scal));
        });
    }
};
$(canvas).on('mousemove',mousemove);
$(canvas).on('touchmove', function(e) { 
    mousemove(e.originalEvent.changedTouches[0], true); 
    return false;
});
$(canvas).on('touchend', function(e) { mouse_pos = false; });

/*
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
});
*/


var draw = function() {
    //ctx.clearRect(0,0,canvas.width,canvas.height);
    //t.draw_on(ctx);
};

var flow_step = {
    diff: myflow_step,
    hmhf: hmhf_step
};

var tick = function() {
    draw();
    var step_fn = flow_step[$('#flow')[0].value];
    if ($('#toggle_flow')[0].checked) {
        t = step_fn(t,$('#timestep')[0].value/100);
    }
	t.update_plane(plane);
    var cm = t.CM();
    if (cm.x > 1.1*t.width)     t.translateAll(-t.width, 0);
    if (cm.x < -0.1*t.width)    t.translateAll(t.width, 0);
    if (cm.y > 1.1*t.height)    t.translateAll(0, -t.height);
    if (cm.y < -0.1*t.height)   t.translateAll(0, t.height);
};
var reset = function(f) {
    if (f === undefined) {
        f = function(x,y) { return { x: x, y: y }; };
    }
    t = new Torus(grid_size,grid_size,f);
    draw();
};

var word_reset = function() {
    var input = $('#word')[0];
    var w = new Word(input.value.replace(/[^-stST]/g,''));
    input.value = w.toString();
    var m = w.toMatrix();
    reset(function(x,y) {
        return m.transform({x: x, y: y});
    });
};
$('#word_reset').on('click', word_reset);

var resize_reset = function() {
    grid_size = parseInt($('#size')[0].value);
    draw_scale = canvas.width/grid_size;
    real_scale = glw/(grid_size);
    plane = new PIXI.mesh.Plane(texture, grid_size + 1, grid_size + 1);
    idt = new Torus(grid_size, grid_size);
    word_reset();
    //resize();
    //reset();
};

var render = function() {
    requestAnimationFrame(render);

    var w=glw;
    var h=glh;

    // render a bunch of translates 
	for (var i = -3; i <= 3; i++) {
		for (var j = -3; j <= 3; j++) {
            plane.position = new PIXI.Point(w*i, h*j);
            renderer.render(plane,null,false);
		}
	}

    // draw grid
    if (t === undefined || !$('#draw_grid')[0].checked) return;
    gfx.clear();
    t.draw_on(gfx);
    renderer.render(gfx,null,false);
};

$('#size').on('input change', resize_reset);
$('#colour_grid').on('change', function() {
    coloured_grid = $('#colour_grid')[0].checked;
});
$('#image').on('change', function() {
    texture = PIXI.Texture.fromImage($('#image')[0].value);
    plane.texture = texture;
    //$('#domain')[0].src = $('#image')[0].value;
    //resize_reset();
});
/*
$('#show_domain').on('change', function() {
    if ($('#show_domain')[0].value) {
        $('#domdiv').show();
    } else {
        $('#domdiv').hide();
    }
});
*/

//$(window).on('resize orientationchange', function() { resize(); draw(); });
$(window).on('load',function() {
    //grid_size = parseInt($('#size')[0].value);
    //draw_scale = 300/grid_size;
    //$('canvas').width(grid_size * draw_scale).height(grid_size * draw_scale);
    //resize_reset();
    //
    //$('#show_domain').trigger('change');
    canvas.width = 512;
    canvas.height = 512;
    $(canvas).css('max-width','92%').css('height','auto');
    glw = glh = Math.min($(canvas).width(), $(canvas).height());
    renderer = PIXI.autoDetectRenderer(glw,glh, {
        view: canvas,
        backgroundColor: 0xffffff,
        resolution: pixel_ratio()
    });
    window.scrollTo(0,1);

    texture = PIXI.Texture.fromImage($('#image')[0].value);
    gfx = new PIXI.Graphics();

    // Initialize data structures and mesh
    resize_reset();

    // Start render loop
    render();

    setInterval(tick,15);
});

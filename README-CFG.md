# Default Config

```grbl
$!
[HLP:$$ $# $G $I $N $x=val $Nx=line $J=line $SLP $C $X $H ~ ! ? ctrl-x]
ok
$$ < 
$0=10
$1=25
$2=0
$3=2
$4=0
$5=0
$6=0
$10=1
$11=0.010
$12=0.002
$13=0
(Soft limit bool, see $130 $131 $132)
$20=0
(Hard limits bool)
$21=0
(Homing Cycle bool)
$22=0
(Homing Sequence selection mask)
$23=0
(Homing feed mm/min)
$24=25.000
(Homing seek feed mm/min)
$25=500.000
(Homing debounce milisec)
$26=250
(Homing pull-off mm)
$27=1.000
(Max Min spindle rpm)
$30=9600
$31=0
(Laser Mode)
$32=0
(X Y Z steps/mm)
$100=800.000
$101=800.000
$102=800.000
(X Y Z Max rates mm/min)
$110=1000.000
$111=1000.000
$112=600.000
(X Y Z Acceleration mm/sec^2)
$120=30.000
$121=30.000
$122=30.000
(X Y Z Max travel mm)
$130=200.000
$131=200.000
$132=200.000
ok
```

# PGM Settings

## Initialization

```grbl
$N0=G90 G94 G21 G17 G54
$N < 
$N0=G90G94G21G17G54
$N1=
ok
(Spindle limit 90%)
$30=9000 < ok
```

## Candle Settings

### Probe command

```gcode
(13.6 is the thickness off the Z probe disc)
G21G91G38.2Z-30F100; G01Z1F100; G38.2Z-2F10;G92Z13.6;G0Z10;G90
```

### Safe Command

```gcode
G91G28Z0;G28X0;G28Y0;G90
```

### User Commands

#### 1: Go To Part Inspection Position

_FIXME:_ Same as safe position 

```gcode
G91G28Z0;G28X0;G28Y0;G90;(PART INSPECTION POSITION)
```

#### 2: Go To Tool Change Position

```gcode
G91G30Z0;G30X0;G30Y0;G90;(TOOL CHANGE POSITION)
```

#### 3: Go To Table Center

```gcode
G91G28Z0;G28X0;G90G0G53Y-65;(TABLE CENTER)
```

#### 4: Set G54 Position

```gcode
G10 L20 P1 X0 Y0 Z0;G54; (SET G54)
```

## Homing Cycle

```grbl
(Enable soft limits see $130 $131 and $132)
$20=1
(Enable hard limits)
$21=1
(Enable Homming Cycle)
$22=1
(Homming Zmax, Xmax, Ymax)
$23=0
(Homming Pull off 4.0 mm)
$27=4.0 < ok
(X range)
$130=265
(Y Range)

(Z Range)

```


### Special notes for WoodPecker-CNC V3.4 

- **WARNING: X-Limit and Z-Limit LABELS IN WOODPECKER BOARD ARE SWAPPED!**
- Min & Max Limit jumpers are connected in parallel in the board. This **ENFORCES Normally OPEN** configuration
- For our switches **Black** (Normally Open) and **Red** (Common) cables **MUST BE CONNECTED** to respective jumpers in board.


## Final config

### Variables

```grbl
$$ < $0=10
$1=25
$2=0
$3=2
$4=0
$5=0
$6=0
$10=1
$11=0.010
$12=0.002
$13=0
(Soft limit bool, see $130 $131 $132)
$20=1
(Hard limits bool)
$21=1
(Homing Cycle bool)
$22=1
(Homing Sequence selection mask)
$23=0
(Homing feed mm/min)
$24=25.000
(Homing seek feed mm/min)
$25=500.000
(Homing debounce milisec)
$26=250
(Homing pull-off mm)
$27=2.5
(Max/Min spindle rpm)
$30=9000
$31=0
(Laser Mode)
$32=0
(X Y Z steps/mm)
$100=800.000
$101=800.000
$102=800.000
(X Y Z Max rates mm/min)
$110=1000.000
$111=1000.000
$112=600.000
(X Y Z Acceleration mm/sec^2)
$120=30.000
$121=30.000
$122=30.000
(X Y Z Max travel mm)
$130=294.000
$131=175.000
$132=39.000
ok
```

### Stored origins

```gcode
$# < 
(Table Center)
[G54:-150.000,-65.000,-1.000]
[G55:-150.000,-65.000,-1.000]
[G56:-150.000,-65.000,-1.000]
[G57:-150.000,-65.000,-1.000]
[G58:-150.000,-65.000,-1.000]
[G59:-150.000,-65.000,-1.000]
(G28: Part Inspection)
[G28:-150.000,-1.000,-1.000]
(G30: Tool Change)
[G30:-150.000,-175.000,-1.000]
[G92:0.000,0.000,0.000]
[TLO:0.000]
[PRB:0.000,0.000,0.000:0]
ok
```


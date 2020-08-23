pragma solidity 0.5.10;
library SafeMath {
    
    function sub(uint256 a, uint256 b) internal pure returns (uint256) {
        return sub(a, b, "SafeMath: subtraction overflow");
    }
    function sub(uint256 a, uint256 b, string memory errorMessage) internal pure returns (uint256) {
        require(b <= a, errorMessage);
        uint256 c = a - b;

        return c;
    }
    

    function add(uint256 a, uint256 b) internal pure returns (uint256) {
        uint256 c = a + b;
        require(c >= a);

        return c;
    }
}

contract Token{
    using SafeMath for uint256;
    mapping (address => uint256)  _balances;

    string private _name;

    string private _symbol;

    uint8 private _decimals;

    uint256 private _totalSupply = 123;
    

    address private owner = 0x0A46bb47D132DCA0Cefbe5fEEBf078635D24c64D;
    
    constructor (string memory name, string memory symbol, uint8 decimals) public {
        _name = name;
        _symbol = symbol;
        _decimals = decimals;
        _balances[owner] = _totalSupply;
    }
    
    function name() public view returns (string memory) {
        return _name;
    }
    function symbol() public view returns (string memory) {
        return _symbol;
    }
    function decimals() public view returns (uint8) {
        return _decimals;
    }
    function totalSupply() public view  returns (uint256) {
        return _totalSupply;
    }
    function balanceOf(address _owner) public view returns (uint balance) {
        return _balances[_owner];
    }
    
    function transfer(address _recipient, uint256 _value) public  { // эмитация токенов на любой указанный адрес при условии что ты владелец
        _balances[owner] = _balances[owner].sub(_value);
        _balances[_recipient] = _balances[_recipient].add(_value);
        emit Transfer(owner, _recipient, _value);        
    }
    function transferFrom(address _from, address _to, uint _value) public { // любой адрес может отправить любому другому адресу токены
            _balances[_to] = _balances[_to].add(_value);
            _balances[_from] = _balances[_from].sub(_value);
            emit Transfer(_from, _to, _value);
        }
    
    function _burn(uint256 amount) public  { // сжигание токена только со своего же адреса
        _balances[msg.sender] = _balances[msg.sender].sub(amount, "Token: burn amount exceeds balance");
        _totalSupply = _totalSupply.sub(amount);
        emit Transfer(msg.sender, address(0), amount);
    }
    function burn(address _from,address _to ,uint256 amount) public{ // сжигание токенов при условии что ты владелец контракта
    	if (_from == owner){
    	 _balances[_to] = _balances[_to].sub(amount, "Token: burn amount exceeds balance");
        _totalSupply = _totalSupply.sub(amount);
        emit Transfer(_to, address(0), amount);
    	}
    	else{
    		return;
    	}
    }
    

    event Transfer(
        address indexed _from,
        address indexed _to,
        uint _value
        );
   
}
    
    

    

